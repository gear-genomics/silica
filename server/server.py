#! /usr/bin/env python

import os
import uuid
import re
import subprocess
import argparse
import json
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
app.config['BASEURL'] = '/silica'
app.static_folder = app.static_url_path = os.path.join(SILICAWS, "../client/src/static")

def allowed_file(filename):
   return '.' in filename and filename.rsplit('.', 1)[1].lower() in set(['fasta', 'fa', 'json', 'csv', 'txt'])

uuid_re = re.compile(r'(^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})-{0,1}([ap]{0,1})([cj]{0,1})$')
def is_valid_uuid(s):
   return uuid_re.match(s) is not None

@app.route('/api/v1/download/<uuid>')
def download(uuid):
   if is_valid_uuid(uuid):
      ma = uuid_re.match(uuid)
      filename = "silica_" + ma.group(1)
      downName = "silica"
      if ma.group(2) == 'a':
         filename += "_amplicons"
         downName += "_amplicons"
      else:
         filename += "_primer"
         downName += "_primer"
      if ma.group(3) == 'c':
         filename += ".csv"
         downName += ".csv"
         miType = 'text/csv'
      else:
         filename += ".json"
         downName += ".json"
         miType = 'text/json'
      if allowed_file(filename):
         sf = os.path.join(app.config['UPLOAD_FOLDER'], uuid[0:2])
         if os.path.exists(sf):
            if os.path.isfile(os.path.join(sf, filename)):
               return send_file(os.path.join(sf, filename), mimetype=miType, as_attachment=True, attachment_filename=downName)
   return "File does not exist!"


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
    outfile = os.path.join(sf, "silica_" + uuidstr + "_amplicons")
    prfile = os.path.join(sf, "silica_" + uuidstr + "_primer")
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

                slexe = os.path.join(app.config['SILICA'], "./src/silica")
                return_code = call([slexe, '-g', genome, '-o', outfile, '-p', prfile,
                                           '--maxProdSize', setAmpSize, '--cutTemp', setTmCutoff,
                                           '--kmer', setKmer, '--distance', setEDis,
                                           '--cutoffPenalty', setCutoffPen, '--penaltyTmDiff', setPenTmDiff,
                                           '--penaltyTmMismatch', setPenTmMismatch, '--penaltyLength', setPenLength,
                                           '--monovalent', setCtmMv, '--divalent', setCtmDv,
                                           '--dna', setCtmDNA, '--dntp', setCtmDNTP,
                                           '-f', 'jsoncsv', ffaname], stdout=log, stderr=err)

    if return_code != 0:
        errInfo = "!"
        with open(errfile, "r") as err:
            errInfo = ": " + err.read()
        return jsonify(errors = [{"title": "Error in running silica" + errInfo}]), 400
    alldata = '{'
    with open(outfile + ".json") as amp:
        ampData = amp.read()
        if len(ampData) > 0:
            alldata += '"amplicon":' + ampData
        else:
            alldata += '"amplicon":[]'
    alldata += ','
    with open(prfile + ".json") as pr:
        prData = pr.read()
        if len(prData) > 0:
            alldata += '"primer":' + prData
        else:
            alldata += '"primer":[]'
    alldata += ', "uuid":"' + uuidstr + '"'
    alldata += '}'
    return jsonify(data = json.loads(alldata)), 200


@app.route('/api/v1/results/<uuid>', methods = ['GET', 'POST'])
def results(uuid):
    if is_valid_uuid(uuid):
        sf = os.path.join(app.config['UPLOAD_FOLDER'], uuid[0:2])
        if os.path.exists(sf):
            ampfilename = "silica_" + uuid + "_amplicons.json";
            if allowed_file(ampfilename):
                if os.path.isfile(os.path.join(sf, ampfilename)):
                    primfilename = "silica_" + uuid + "_primer.json";
                    if allowed_file(primfilename):
                        if os.path.isfile(os.path.join(sf, primfilename)):
                            alldata = '{'
                            with open(os.path.join(sf, ampfilename)) as amp:
                                ampData = amp.read()
                                if len(ampData) > 0:
                                    alldata += '"amplicon":' + ampData
                                else:
                                    alldata += '"amplicon":[]'
                            alldata += ','
                            with open(os.path.join(sf, primfilename)) as pr:
                                prData = pr.read()
                                if len(prData) > 0:
                                    alldata += '"primer":' + prData
                                else:
                                    alldata += '"primer":[]'
                            alldata += ', "uuid":"' + uuid + '"'
                            alldata += '}'
                            return jsonify(data = json.loads(alldata)), 200
    return jsonify(errors = [{"title": "Link outdated or invalid!"}]), 400

@app.route('/api/v1/genomeindex', methods=['POST'])
def genomeind():
    return send_from_directory(os.path.join(SILICAWS, "../fm"),"genomeindexindex.json"), 200

def onlyFloat(txt):
    onlyNumbDC = re.compile('[^0-9,.\-]')
    txt = onlyNumbDC.sub( '', txt)
    txt = txt.replace(',', '.')
    return txt

def onlyInt(txt):
    onlyNumb = re.compile('[^0-9\-]')
    return onlyNumb.sub( '', txt)

if __name__ == '__main__':
    app.run(host = '0.0.0.0', port=3300, debug = True, threaded=True)

                                                                        
