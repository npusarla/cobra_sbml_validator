from gzip import GzipFile
from bz2 import decompress as bz2_decompress
from tempfile import NamedTemporaryFile
from json import dumps
from os import mkdir, unlink, path
from os.path import isdir, isfile, join
from Bio import Entrez, SeqIO
import argparse 
import re
from warnings import catch_warnings
from codecs import getreader
import time

from six import BytesIO, StringIO, iteritems
import jsonschema

import cobra
from cobra.core.gene import parse_gpr
from cobra.manipulation import check_mass_balance, check_reaction_bounds, \
    check_metabolite_compartment_formula

from libsbml import SBMLValidator

from flask import Flask, session, render_template, request, jsonify
from celery import Celery
from flask_session import Session
app = Flask(__name__, template_folder='.')
SESSION_TYPE = 'redis'
app.config.from_object(__name__)
Session(app)

validator_form = path.join(path.abspath(path.dirname(__file__)),
                           "validator_form.html")

app = Flask(__name__)
app.config['CELERY_BROKER_URL'] = 'redis://localhost:6379/0'
app.config['CELERY_RESULT_BACKEND'] = 'redis://localhost:6379/0'

#initialize celery
celery = Celery(app.name, broker=app.config['CELERY_BROKER_URL'])
celery.conf.update(app.config)

FOLDER_NAME = 'genbank'
if not isdir(FOLDER_NAME):
    mkdir(FOLDER_NAME)

@celery.task
def load_JSON(contents):
    """returns model, [model_errors], "parse_errors" or None """
    errors = []
    try:
        model_json = cobra.io.json.json.load(getreader("utf-8")(contents))
    except ValueError as e:
        return None, errors, "Invalid JSON: " + str(e)
    try:
        model = cobra.io.json._from_dict(model_json)
    except Exception as e:
        errors.append("Invalid model: " + str(e))
        model = None
    try:
        jsonschema.validate(model_json, cobra.io.json.json_schema)
    except jsonschema.ValidationError as e:
        # render an infomrative error message
        if len(e.absolute_path) > 0:
            error_msg = "Error in "
            for i in e.absolute_path:
                if isinstance(i, int):
                    error_msg = error_msg.rstrip(".") + "[%d]." % i
                else:
                    error_msg += str(i) + "."
            errors.append(error_msg.rstrip(".") + ": " + e.message)
        else:
            errors.append(e.message)
    return model, errors, None

@celery.task
def load_SBML(contents, filename):
    """returns model, [model_errors], "parse_errors" or None """
    try:  # this function fails if a model can not be created
        model, errors = cobra.io.sbml3.validate_sbml_model(
            contents, check_model=False)  # checks are run later
    except cobra.io.sbml3.CobraSBMLError as e:
        return None, [], str(e)
    else:
        return model, errors, None


@celery.task
def run_libsbml_validation(contents, filename):
    if filename.endswith(".gz"):
        filename = filename[:-3]
    elif filename.endswith(".bz2"):
        filename = filename[:-4]
    with NamedTemporaryFile(suffix=filename, delete=False) as outfile:
        outfile.write(contents.read())
        contents.seek(0)  # so the buffer can be re-read
    validator = SBMLValidator()
    validator.validate(str(outfile.name))
    unlink(outfile.name)
    errors = []
    for i in range(validator.getNumFailures()):
        failure = validator.getFailure(i)
        if failure.isWarning():
            continue
        errors.append("L%d C%d: %s" % (failure.getLine(), failure.getColumn(),
                                       failure.getMessage()))
    return errors

@celery.task
def decompress_file(body, filename):
    """returns BytesIO of decompressed file"""
    if filename.endswith(".gz"):
        # contents = zlib.decompress(body, 16 + zlib.MAX_WBITS)
        zip_contents = BytesIO(body)
        with GzipFile(fileobj=zip_contents, mode='rb') as zip_read:
            try:
                contents = BytesIO(zip_read.read())
            except (IOError, OSError) as e:
                return None, "Error decompressing gzip file: " + str(e)
        zip_contents.close()
    elif filename.endswith(".bz2"):
        try:
            contents = BytesIO((bz2_decompress(body)))
        except IOError as e:
            return None, "Error decompressing bz2 file: " + str(e)
    else:
        contents = BytesIO((body))
    return contents, None


@celery.task
def validate_model(model):
    errors = []
    warnings = []
    errors.extend(check_reaction_bounds(model))
    errors.extend(check_metabolite_compartment_formula(model))
    # test gpr
    for reaction in model.reactions:
        try:
            parse_gpr(reaction.gene_reaction_rule)
        except SyntaxError:
            errors.append("reaction '%s' has invalid gpr '%s'" %
                          (reaction.id, reaction.gene_reaction_rule))
    # test mass balance
    for reaction, balance in iteritems(check_mass_balance(model)):
        # check if it's a demand or exchange reaction
        if len(reaction.metabolites) == 1:
            warnings.append("reaction '%s' is not balanced. Should it "
                            "be annotated as a demand or exchange "
                            "reaction?" % reaction.id)
        elif "biomass" in reaction.id.lower():
            warnings.append("reaction '%s' is not balanced. Should it "
                            "be annotated as a biomass reaction?" %
                            reaction.id)
        else:
            warnings.append("reaction '%s' is not balanced for %s" %
                            (reaction.id, ", ".join(sorted(balance))))

    # try solving
    solution = model.optimize()
    if solution.status != "optimal":
        errors.append("model can not be solved (status '%s')" %
                      solution.status)
        return {"errors": errors, "warnings": warnings}

    # if there is no objective, then we know why the objective was low
    if len(model.objective) == 0:
        warnings.append("model has no objective function")
    elif solution.f <= 0:
        warnings.append("model can not produce nonzero biomass")
    elif solution.f <= 1e-3:
        warnings.append("biomass flux %s too low" % str(solution.f))
    if len(model.objective) > 1:
        warnings.append("model should only have one reaction as the objective")

    return {"errors": errors, "warnings": warnings, "objective": solution.f}


def gen_filepath(accession):
    return join(FOLDER_NAME, accession + '.gb')

# def write_error(self, status_code, reason="", **kwargs):
#     self.write(reason)


@app.route('/status/<task_id>')
def taskstatus(task_id):
    task = handle_uploaded_file.AsyncResult(task_id)
    if task.state == 'PENDING':
        # job did not start yet
        response = {
            'state': task.state,
            'status': 'Pending...'
        }
    elif task.state != 'FAILURE':
        response = {
            'state': task.state,
            'model_warnings': task.info.get('model_warnings', []),
            'model_errors': task.info.get('model_errors', []),
            'status': task.info.get('status', '')
        }
    else:
        # something went wrong in the background job
        response = {
            'state': task.state,
            'status': str(task.info),  # this is the exception raised
        }
    return jsonify(response)


@app.route('/upload', methods=['POST'])
def upload():

    print(request.files)
    fileinfo = request.files['file']
    filename = fileinfo.filename

    #data = request.form["ncbi_accession"]

    # option 1: synchronous, < 100ms
    # t = time.time()
    # contents, error = decompress_file(fileinfo.body, filename)
    # print(f'took {time.time() - t}')

    # option 2: synchronous with celery, < 100ms
    # t = time.time()
    # result = decompress_file.apply_async(args=(fileinfo.body, filename))
    # contents, error = result.get()
    # print(f'took {time.time() - t}')

    # option 2: async within a HTTP query, < 2 seconds
    # contents, error = yield executor.submit(
    #      decompress_file_fn, fileinfo["body"], filename)

    # option 3: task in background (actually synchronous), > 2 seconds
    task = handle_uploaded_file.apply_async(args=(fileinfo.body, filename))
    # TODO return task.id to user
    return jsonify({'task_id': task.id})



    # return jsonify({'errors': ['test error'], 'warnings': ['test warning']})

    # if error:
        #error handling?

     # syntax validation
     # if the model can't be loaded from the file it's considered invalid

     # if not explicitly JSON, assumed to be SBML
     # warnings = []
     # if filename.endswith(".json") or filename.endswith(".json.gz") or \
     #        filename.endswith(".json.bz2"):
     #    model, errors, parse_errors = \
    #         yield executor.submit(load_JSON, contents)

    # else:
    #     model, errors, parse_errors = \
    #         yield executor.submit(load_SBML, contents, filename)
    #     libsbml_errors = yield executor.submit(
    #         run_libsbml_validation, contents, filename)
    #     warnings.extend("(from libSBML) " + i for i in libsbml_errors)

    # # if parsing failed, then send the error
    # if parse_errors:
    #     self.send_error(415, reason=parse_errors)
    #     return
    # elif model is None:  # parsed, but still could not generate model
    #     self.finish({"errors": errors, "warnings": warnings})
    #     return

    # # model validation
    # result = yield executor.submit(validate_model, model)
    # result["errors"].extend(errors)
    # result["warnings"].extend(warnings)
    
    # #import IPython; IPython.embed()
    
    # gb_filepath = gen_filepath(data)
    # if not isfile(gb_filepath):
    #     dl = Entrez.efetch(db='nuccore', id=data, rettype='gbwithparts',
    #                     retmode='text')
    #     with open(gb_filepath, 'w') as outfile:
    #         outfile.write(dl.read())
    #     dl.close()
    #     print('------------------ DONE writing') 
    

    # #pseudocode
    # gb_seq = SeqIO.read(gb_filepath, 'genbank') 
    # locus_list = []
    # for feature in gb_seq.features:
    #    if feature.type == 'CDS':
    #        locus_list.append(feature.qualifiers['locus_tag'][0])

    # #backup list 
    # gene_list = []
    # for feature in gb_seq.features:
    #     if feature.type == 'CDS':
    #         for i in feature:
    #             if i is 'gene':
    #                 gene_list.append(i[0])

    # #pseudocode for checking the genes
    # model_genes = model.genes

    # #import IPython; IPython.embed()
           
    # model_genes_id = []
    # for gene in model_genes:
    #     model_genes_id.append(gene.id)

    # badGenes = list(set(model_genes_id) - set(locus_list))
    # #backup search 
    # badGenes2 = list(set(badGenes) - set(gene_list))
    # genesToChange = list(set(badGenes).intersection(gene_list))
    # #overallList = list(set(badGenes).intersection(badGenes2))
    # result['errors'].extend([x + ' is not a valid gene' for x in badGenes2])

    # if len(genesToChange) != 0:
    #     result['warnings'].extend(['Change ' + x + ' to locus tag name' for x in genesToChange])


    # self.finish(result)
    
   
    # print('------------------ DONE getting genes')
    
    # #self.finish({ 'status': 'downloaded' })
    
@app.route('/')
def get():
    return render_template('validator_form.html')

if __name__ == "__main__":
    app.run(threaded=True, debug=True)

    # import argparse
    # parser = argparse.ArgumentParser(
    #     description="web-based validator for COBRA models in SBML and JSON")
    # parser.add_argument("--port", type=int, default=5000)
    # parser.add_argument("--prefix", default="")
    # parser.add_argument("--debug", action="store_true")

    # args = parser.parse_args()
    
    # prefix = args.prefix
    # if len(prefix) > 0 and not prefix.startswith("/"):
    #     prefix = "/" + prefix

    
    # run_standalone_server(
    #     prefix=prefix,
    #     port=args.port,
    #     debug=args.debug)
    
