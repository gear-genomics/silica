import os
import uuid
import re
import subprocess
import argparse
from subprocess import call
from flask import Flask, send_file, flash, render_template, request, redirect, url_for
from werkzeug.utils import secure_filename


SILICAWS = os.path.dirname(os.path.abspath(__file__))

app = Flask(__name__)
app.config['SILICA'] = os.path.join(SILICAWS, "..")
app.config['UPLOAD_FOLDER'] = os.path.join(app.config['SILICA'], "data")
app.config['MAX_CONTENT_LENGTH'] = 8 * 1024 * 1024   #maximum of 8MB
app.config['BASEURL'] = '/silica'


def allowed_file(filename):
   return '.' in filename and filename.rsplit('.', 1)[1].lower() in set(['fasta', 'fa', 'txt'])

uuid_re = re.compile(r'^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$')
def is_valid_uuid(s):
   return uuid_re.match(s) is not None

@app.route('/download/<uuid>')
def download(uuid):
   if is_valid_uuid(uuid):
      filename = "silica_" + uuid + "_amplicons.txt";
      if allowed_file(filename):
         sf = os.path.join(app.config['UPLOAD_FOLDER'], uuid[0:2])
         if os.path.exists(sf):
            if os.path.isfile(os.path.join(sf, filename)):
               return send_file(os.path.join(sf, filename), attachment_filename=filename)
   return "File does not exist!"

@app.route('/upload', methods = ['GET', 'POST'])
def upload_file():
   if request.method == 'POST':
      uuidstr = str(uuid.uuid4())

      # Get subfolder
      sf = os.path.join(app.config['UPLOAD_FOLDER'], uuidstr[0:2])
      if not os.path.exists(sf):
         os.makedirs(sf)

      # Load Test Data
      if request.form['submit'] == 'Load Example Data':
         testData = '>FGA_f\nGCCCCATAGGTTTTGAACTCA\n>FGA_r\nTGATTTGTCTGTAATTGCCAGC'
         return render_template('upload.html', baseurl = app.config['BASEURL'], example = testData)

      # Fasta file
      primerData = request.form['fastaText']
      primerFile = 'direct.fa'
      if 'fasta' in request.files:
         ffa = request.files['fasta']
         if ffa.filename != '':
            if allowed_file(ffa.filename):
               primerData += '\n' + ffa.read()
               primerFile = secure_filename(ffa.filename)
      if  primerData == '':
         error = "Please provide a set of primers!"
         return render_template('upload.html', baseurl = app.config['BASEURL'], error = error)
      ffaname = os.path.join(sf, "silica_" + uuidstr + "_" + primerFile)
      primerData = primerData.replace('\r\n','\n')
      with open(ffaname, "w") as primFile:
         primFile.write(primerData)
    #  ffa.save(ffaname)

      # Genome
      val = request.form.get("submit", "None provided")
      genome = request.form['genome']
      if genome == '':
         error = "Please select a genome!"
         return render_template('upload.html', baseurl = app.config['BASEURL'], error = error)
      genome = os.path.join(app.config['SILICA'], "fm", genome)

      # Run Rscript
      outfile = os.path.join(sf, "silica_" + uuidstr + "_amplicons.txt")
      prfile = os.path.join(sf, "silica_" + uuidstr + "_primer.txt")
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
                  error = 'Mimnimal Primer Tm must be >= 30&deg;C'
                  return render_template('upload.html', baseurl = app.config['BASEURL'], error = error)
               setKmer = onlyInt(request.form['setKmer'])
               param.write("kmer=" + setKmer + '\n')
               if int(setKmer) < 15:
                  error = 'Number of bp Used to Seach for Matches must be > 14 bp'
                  return render_template('upload.html', baseurl = app.config['BASEURL'], error = error)
               setEDis = onlyInt(request.form['setDist'])
               param.write("distance=" + setEDis + '\n')
               if int(setEDis) != 0 and int(setEDis) != 1:
                  error = 'Maximal Allowed Number of Mutations must be 0 or 1'
                  return render_template('upload.html', baseurl = app.config['BASEURL'], error = error)
               setCutoffPen = onlyFloat(request.form['setCutoffPen'])
               param.write("cutoffPenalty=" + setCutoffPen + '\n')
               if float(setCutoffPen) < 0.0 and int(setCutoffPen) != -1:
                  error = 'Keep Only PCR Products with Penalty Below must be > 0.0 or -1'
                  return render_template('upload.html', baseurl = app.config['BASEURL'], error = error)
               setPenTmDiff = onlyFloat(request.form['setPenTmDiff'])
               param.write("penaltyTmDiff=" + setPenTmDiff + '\n')
               if float(setPenTmDiff) < 0.0:
                  error = 'Penalty Factor for Single Primer Tm Mismatch must be >= 0.0'
                  return render_template('upload.html', baseurl = app.config['BASEURL'], error = error)
               setPenTmMismatch = onlyFloat(request.form['setPenTmMismatch'])
               param.write("penaltyTmMismatch=" + setPenTmMismatch + '\n')
               if float(setPenTmMismatch) < 0.0:
                  error = 'Penalty Factor for Tm Mismatch of Primers in a Pair must be >= 0.0'
                  return render_template('upload.html', baseurl = app.config['BASEURL'], error = error)
               setPenLength = onlyFloat(request.form['setPenLength'])
               param.write("penaltyLength=" + setPenLength + '\n')
               if float(setPenLength) < 0.0:
                  error = 'Penalty Factor for PCR Product Length must be >= 0.0'
                  return render_template('upload.html', baseurl = app.config['BASEURL'], error = error)
               setCtmMv = onlyFloat(request.form['setCtmMv'])
               param.write("monovalent=" + setCtmMv + '\n')
               if float(setCtmMv) < 0.0:
                  error = 'Concentration of Monovalent Ions must be >= 0.0 mMol'
                  return render_template('upload.html', baseurl = app.config['BASEURL'], error = error)
               setCtmDv = onlyFloat(request.form['setCtmDv'])
               param.write("divalent=" + setCtmDv + '\n')
               if float(setCtmDv) < 0.0:
                  error = 'Concentration of Divalent Ions must be >= 0.0 mMol'
                  return render_template('upload.html', baseurl = app.config['BASEURL'], error = error)
               setCtmDNA = onlyFloat(request.form['setCtmDNA'])
               param.write("dna=" + setCtmDNA + '\n')
               if float(setCtmDNA) < 0.1:
                  error = 'Concentration of Annealing(!) Oligos must be >= 0.1 nMol'
                  return render_template('upload.html', baseurl = app.config['BASEURL'], error = error)
               setCtmDNTP = onlyFloat(request.form['setCtmDNTP'])
               param.write("dntp=" + setCtmDNTP + '\n')
               if float(setCtmDNTP) < 0.0:
                  error = 'Concentration of the Sum of All dNTPs must be >= 0.0 mMol'
                  return render_template('upload.html', baseurl = app.config['BASEURL'], error = error)

               slexe = os.path.join(app.config['SILICA'], "./src/silica")
               return_code = call([slexe, '-g', genome, '-o', outfile, '-p', prfile,
                                          '--maxProdSize', setAmpSize, '--cutTemp', setTmCutoff,
                                          '--kmer', setKmer, '--distance', setEDis,
                                          '--cutoffPenalty', setCutoffPen, '--penaltyTmDiff', setPenTmDiff,
                                          '--penaltyTmMismatch', setPenTmMismatch, '--penaltyLength', setPenLength,
                                          '--monovalent', setCtmMv, '--divalent', setCtmDv,
                                          '--dna', setCtmDNA, '--dntp', setCtmDNTP,
                                          ffaname], stdout=log, stderr=err)
      if return_code != 0:
         error = "Error in running Silica!"
         return render_template('upload.html', baseurl = app.config['BASEURL'], error = error)

      # Send download pdf
      return redirect(app.config['BASEURL'] + "/download/" + uuidstr, code=302)
   return render_template('upload.html', baseurl = app.config['BASEURL'])

@app.route("/")
def submit():
    return render_template("upload.html", baseurl = app.config['BASEURL'])

def onlyFloat(txt):
    onlyNumbDC = re.compile('[^0-9,.\-]')
    txt = onlyNumbDC.sub( '', txt)
    txt = txt.replace(',', '.')
    return txt

def onlyInt(txt):
    onlyNumb = re.compile('[^0-9\-]')
    return onlyNumb.sub( '', txt)

if __name__ == '__main__':
   parser = argparse.ArgumentParser(description='Silica App')
   parser.add_argument('-b', '--baseurl', type=str, required=False, default="", metavar="", dest='baseurl', help='baseurl')
   parser.add_argument('-d', '--debug', required=False, default=False, action='store_true', dest='debug', help='debug mode')
   parser.add_argument('-t', '--threaded', required=False, default=False, action='store_true', dest='threaded', help='threaded')
   parser.add_argument('-s', '--host', type=str, required=False, default="0.0.0.0", metavar="0.0.0.0", dest='host', help='host')
   parser.add_argument('-p', '--port', type=int, required=False, default=3300, metavar="3300", dest='port', help='port')
   args = parser.parse_args()
   app.config['BASEURL'] = args.baseurl
   app.run(host = args.host, port = args.port, debug = args.debug, threaded= args.threaded)
