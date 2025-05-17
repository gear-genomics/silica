#! /usr/bin/env python

import os
import errno
import uuid
import re
import subprocess
import argparse
import json
import gzip
from subprocess import call
from flask import Flask, send_file, flash, send_from_directory, request, redirect, url_for, jsonify
from flask_cors import CORS
from werkzeug.utils import secure_filename

app = Flask(__name__)
CORS(app)
SILICAWS = os.path.dirname(os.path.abspath(__file__))

app.config['SILICA'] = os.path.join(SILICAWS, "..")
app.config['UPLOAD_FOLDER'] = os.path.join(app.config['SILICA'], "data")
app.config['MAX_CONTENT_LENGTH'] = 8 * 1024 * 1024   #maximum of 8MB


def allowed_file(filename):
   return '.' in filename and filename.rsplit('.', 1)[1].lower() in set(['fasta', 'fa', 'json', 'csv', 'txt'])


uuid_re = re.compile(r'(^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})-{0,1}([ap]{0,1})([cj]{0,1})$')
def is_valid_uuid(s):
   return uuid_re.match(s) is not None


@app.route('/api/v1/upload', methods=['POST'])
def generate():
    uuidstr = str(uuid.uuid4())

    # Get subfolder
    sf = os.path.join(app.config['UPLOAD_FOLDER'], uuidstr[0:2])
    if not os.path.exists(sf):
        os.makedirs(sf)

    # Fasta file
    primerData = request.form['fastaText']
    primerData = primerData.replace('\r\n','\n')
    if primerData == '':
        return jsonify(errors = [{"title": "Please provide a set of primers!"}]), 400
    ffaname = os.path.join(sf, "silica_" + uuidstr + "_fasta.fa")
    with open(ffaname, "w") as primFile:
        primFile.write(primerData)

    # Genome
    val = request.form.get("submit", "None provided")
    genome = request.form['genome']
    if genome == '':
        return jsonify(errors = [{"title": "Please select a genome!"}]), 400
    genome = os.path.join(app.config['SILICA'], "fm", genome)

    # Run silica
    outfile = os.path.join(sf, "silica_" + uuidstr + ".json.gz")
    paramfile = os.path.join(sf, "silica_" + uuidstr + "_parameter.txt")
    logfile = os.path.join(sf, "silica_" + uuidstr + ".log")
    errfile = os.path.join(sf, "silica_" + uuidstr + ".err")
    with open(logfile, "w") as log:
        with open(errfile, "w") as err:
            with open(paramfile, "w") as param:
                param.write("genome=" + genome + '\n')
                setAmpSize = onlyInt(request.form['setAmpSize'])
                param.write("maxProdSize=" + setAmpSize + '\n')
                setTmCutoff = onlyFloat(request.form['setTmCutoff'])
                param.write("cutTemp=" + setTmCutoff + '\n')
                if float(setTmCutoff) < 30.0:
                    return jsonify(errors = [{"title": "Mimnimal Primer Tm must be >= 30&deg;C"}]), 400
                setKmer = onlyInt(request.form['setKmer'])
                param.write("kmer=" + setKmer + '\n')
                if int(setKmer) < 15:
                    return jsonify(errors = [{"title": "Number of bp Used to Seach for Matches must be > 14 bp"}]), 400
                setEDis = onlyInt(request.form['setDist'])
                param.write("distance=" + setEDis + '\n')
                if int(setEDis) != 0 and int(setEDis) != 1:
                    return jsonify(errors = [{"title": "Maximal Allowed Number of Mutations must be 0 or 1"}]), 400
                setCutoffPen = onlyFloat(request.form['setCutoffPen'])
                param.write("cutoffPenalty=" + setCutoffPen + '\n')
                if float(setCutoffPen) < 0.0 and int(setCutoffPen) != -1:
                    return jsonify(errors = [{"title": "Keep Only PCR Products with Penalty Below must be > 0.0 or -1"}]), 400
                setPenTmDiff = onlyFloat(request.form['setPenTmDiff'])
                param.write("penaltyTmDiff=" + setPenTmDiff + '\n')
                if float(setPenTmDiff) < 0.0:
                    return jsonify(errors = [{"title": "Penalty Factor for Single Primer Tm Mismatch must be >= 0.0"}]), 400
                setPenTmMismatch = onlyFloat(request.form['setPenTmMismatch'])
                param.write("penaltyTmMismatch=" + setPenTmMismatch + '\n')
                if float(setPenTmMismatch) < 0.0:
                    return jsonify(errors = [{"title": "Penalty Factor for Tm Mismatch of Primers in a Pair must be >= 0.0"}]), 400
                setPenLength = onlyFloat(request.form['setPenLength'])
                param.write("penaltyLength=" + setPenLength + '\n')
                if float(setPenLength) < 0.0:
                    return jsonify(errors = [{"title": "Penalty Factor for PCR Product Length must be >= 0.0"}]), 400
                setCtmMv = onlyFloat(request.form['setCtmMv'])
                param.write("monovalent=" + setCtmMv + '\n')
                if float(setCtmMv) < 0.0:
                    return jsonify(errors = [{"title": "Concentration of Monovalent Ions must be >= 0.0 mMol"}]), 400
                setCtmDv = onlyFloat(request.form['setCtmDv'])
                param.write("divalent=" + setCtmDv + '\n')
                if float(setCtmDv) < 0.0:
                    return jsonify(errors = [{"title": "Concentration of Divalent Ions must be >= 0.0 mMol"}]), 400
                setCtmDNA = onlyFloat(request.form['setCtmDNA'])
                param.write("dna=" + setCtmDNA + '\n')
                if float(setCtmDNA) < 0.1:
                    return jsonify(errors = [{"title": "Concentration of Annealing(!) Oligos must be >= 0.1 nMol"}]), 400
                setCtmDNTP = onlyFloat(request.form['setCtmDNTP'])
                param.write("dntp=" + setCtmDNTP + '\n')
                if float(setCtmDNTP) < 0.0:
                    return jsonify(errors = [{"title": "Concentration of the Sum of All dNTPs must be >= 0.0 mMol"}]), 400

                try:
                    return_code = call(['dicey', 'search', '-g', genome, '-o', outfile, '-i', os.path.join(SILICAWS, "../primer3_config/"),
                                               '--maxProdSize', setAmpSize, '--cutTemp', setTmCutoff,
                                               '--kmer', setKmer, '--distance', setEDis,
                                               '--cutoffPenalty', setCutoffPen, '--penaltyTmDiff', setPenTmDiff,
                                               '--penaltyTmMismatch', setPenTmMismatch, '--penaltyLength', setPenLength,
                                               '--monovalent', setCtmMv, '--divalent', setCtmDv,
                                               '--dna', setCtmDNA, '--dntp', setCtmDNTP,
                                               ffaname], stdout=log, stderr=err)
                except OSError as e:
                    if e.errno == errno.ENOENT:
                        return jsonify(errors = [{"title": "Binary dicey not found!"}]), 400
                    else:
                        return jsonify(errors = [{"title": "OSError " + str(e.errno) + " running binary dicey!"}]), 400
    result = gzip.open(outfile).read()
    if result is None:
        datajs = []
        datajs["errors"] = []
    else:
        datajs = json.loads(result)
    datajs['uuid'] = uuidstr
    with open(errfile, "r") as err:
        errInfo = ": " + err.read()
        if len(errInfo) > 3 or return_code != 0:
            if len(errInfo) > 3:
                datajs["errors"] = [{"title": "Error in running silica" + errInfo}] + datajs["errors"]
            if return_code != 0:
                datajs["errors"] = [{"title": "Run Error - Dicey did not return 0"}] + datajs["errors"]
            return jsonify(datajs), 400
    return jsonify(datajs), 200


@app.route('/api/v1/results/<uuid>', methods = ['GET', 'POST'])
def results(uuid):
    if is_valid_uuid(uuid):
        sf = os.path.join(app.config['UPLOAD_FOLDER'], uuid[0:2])
        if os.path.exists(sf):
            sjsfilename = "silica_" + uuid + ".json.gz"
            if os.path.isfile(os.path.join(sf, sjsfilename)):
                result = gzip.open(os.path.join(sf, sjsfilename)).read()
                if result is None:
                    datajs = []
                    datajs["errors"] = []
                else:
                    datajs = json.loads(result)
                datajs['uuid'] = uuid
                with open(os.path.join(sf, "silica_" + uuid + ".err"), "r") as err:
                    errInfo = ": " + err.read()
                    if len(errInfo) > 3:
                        datajs["errors"] = [{"title": "Error in running silica" + errInfo}] + datajs["errors"]
                        return jsonify(datajs), 400
                return jsonify(datajs), 200
    return jsonify(errors = [{"title": "Link outdated or invalid!"}]), 400


@app.route('/api/v1/genomeindex', methods=['POST'])
def genomeind():
    return send_from_directory(os.path.join(SILICAWS, "../fm"),"genomeindexindex.json"), 200


@app.route('/api/v1/health', methods=['GET'])
def health():
    return jsonify(status="OK")


def onlyFloat(txt):
    onlyNumbDC = re.compile(r'[^0-9,.\-]')
    txt = onlyNumbDC.sub( '', txt)
    txt = txt.replace(',', '.')
    return txt


def onlyInt(txt):
    onlyNumb = re.compile(r'[^0-9\-]')
    return onlyNumb.sub( '', txt)


if __name__ == '__main__':
    app.run(host = '0.0.0.0', port=3300, debug = True, threaded=True)
