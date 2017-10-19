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
      filename = "silica_" + uuid + ".txt";
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
      if request.form['submit'] == 'Load Test Data':
         testData = '>FGA_f\nGCCCCATAGGTTTTGAACTCA\n>FGA_r\nTGATTTGTCTGTAATTGCCAGC'
         return render_template('upload.html', baseurl = app.config['BASEURL'], example = testData)


      # Fasta file



      if 'fasta' not in request.files:
         error = "Fasta file missing!"
         return render_template('upload.html', baseurl = app.config['BASEURL'], error = error)
      ffa = request.files['fasta']
      if ffa.filename == '':
         error = "Fasta file missing!"
         return render_template('upload.html', baseurl = app.config['BASEURL'], error = error)
      if not allowed_file(ffa.filename):
         error = "Fasta file has incorrect file type!"
         return render_template('upload.html', baseurl = app.config['BASEURL'], error = error)
      ffaname = os.path.join(sf, "silica_" + uuidstr + "_" + secure_filename(ffa.filename))
      ffa.save(ffaname)

      # Genome
      val = request.form.get("submit", "None provided")
      genome = request.form['genome']
      if genome == '':
         error = "Genome index is missing!"
         return render_template('upload.html', baseurl = app.config['BASEURL'], error = error)
      genome = os.path.join(app.config['SILICA'], "fm", genome)

      # Run Rscript
      outfile = os.path.join(sf, "silica_" + uuidstr + ".txt")
      prfile = os.path.join(sf, "silica_" + uuidstr + ".primer")
      logfile = os.path.join(sf, "silica_" + uuidstr + ".log")
      errfile = os.path.join(sf, "silica_" + uuidstr + ".err")
      with open(logfile, "w") as log:
         with open(errfile, "w") as err:
            slexe = os.path.join(app.config['SILICA'], "./src/silica")
            return_code = call([slexe, '-g', genome, '-o', outfile, '-p', prfile, ffaname], stdout=log, stderr=err)
      if return_code != 0:
         error = "Error in running Silica!"
         return render_template('upload.html', baseurl = app.config['BASEURL'], error = error)

      # Send download pdf
      return redirect(app.config['BASEURL'] + "/download/" + uuidstr, code=302)
   return render_template('upload.html', baseurl = app.config['BASEURL'])

@app.route("/")
def submit():
    return render_template("upload.html", baseurl = app.config['BASEURL'])

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
