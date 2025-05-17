const API_URL = process.env.API_URL
const API_LINK = process.env.API_LINK

var res = ""
var prStart = 0
var amStart = 0
var step = 30

window.data = ""
genomeConv = {"Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz": "",
              "Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz": "ce11",
              "Danio_rerio.GRCz10.dna.toplevel.fa.gz": "danRer10",
              "Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz": "dm6",
              "Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz": "hg19",
              "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz": "hg38",
              "Mus_musculus.GRCm38.dna.primary_assembly.fa.gz": "mm10",
              "Oryzias_latipes.MEDAKA1.dna.toplevel.fa.gz": "oryLat2",
              "Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz": "sacCer3",
              "SARS-CoV-2.NC_045512.2.dna.fa.gz": "wuhCor1"}
geneBro = "https://genome-euro.ucsc.edu/cgi-bin/hgTracks?"

var fileLoad = document.getElementById('fasta');
fileLoad.addEventListener('change', loadFasta, false);
var submitButton = document.getElementById('btn-submit')
submitButton.addEventListener('click', run)
var sampleButton = document.getElementById('btn-example')
sampleButton.addEventListener('click', sampleData)
var helpButton = document.getElementById('btn-help')
helpButton.addEventListener('click', goToHelp)

const saveButton = document.getElementById('btn-save-Json')
saveButton.addEventListener('click', saveJsonFile)
const loadJFile = document.getElementById('inputJsonFile')
loadJFile.addEventListener('change', loadJsonFile, false);

const resultLink = document.getElementById('link-results')
const helpLink = document.getElementById('link-help')
const resultInfo = document.getElementById('result-info')
const resultError = document.getElementById('result-error')
const resultTabs = document.getElementById('result-tabs')
const primerLink = document.getElementById('link-primers')
const targetGenomes = document.getElementById('target-genome')
var selectedGenome = ""
var sectionResults = document.getElementById('results')
var sectionAmpliconLink = document.getElementById('res-amplicons-link')
var sectionAmpliconData = document.getElementById('res-amplicons-data')
var sectionPrimerLink = document.getElementById('res-primers-link')
var sectionPrimerData = document.getElementById('res-primers-data')

var btnAmpUp = document.getElementById('res-amplicons-Up')
btnAmpUp.addEventListener('click', ampUp)
var btnAmpDown = document.getElementById('res-amplicons-Down')
btnAmpDown.addEventListener('click', ampDown)
var btnPrimerUp = document.getElementById('res-primers-Up')
btnPrimerUp.addEventListener('click', primerUp)
var btnPrimerDown = document.getElementById('res-primers-Down')
btnPrimerDown.addEventListener('click', primerDown)

$('#mainTab a').on('click', function (e) {
  e.preventDefault()
  $(this).tab('show')
})

$('#resTab a').on('click', function (e) {
  e.preventDefault()
  $(this).tab('show')
})

var spinnerHtml = '<p>Analysis takes a couple of seconds to run, please be patient.</p><i class="fa fa-spinner fa-pulse fa-3x fa-fw"></i><br /><br />'

document.addEventListener("DOMContentLoaded", function() {
    checkForUUID();
});

function goToHelp() {
    helpLink.click();
}

function loadFasta(f) {
    var file = f.target.files[0];
    if (file) { // && file.type.match("text/*")) {
        var reader = new FileReader();
        reader.onload = function(event) {
            var txt = event.target.result;
            var regEx1 = /\r\n/g;
            txt = txt.replace(regEx1, "\n");
            document.getElementById('fastaText').value = txt
        }
        reader.readAsText(file);
    } else {
        alert("Error opening file");
    }
}

function sampleData () {
    var fasta = document.getElementById('fastaText')
    fasta.value = '>FGA_f\nGCCCCATAGGTTTTGAACTCA\n>FGA_r\nTGATTTGTCTGTAATTGCCAGC\n';
    var selectbox = document.getElementById('genome-select')
    for (var i = 0 ; i < selectbox.options.length ; i++) {
        if (selectbox.options[i].value == 'Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz') {
            selectbox.selectedIndex = i;
        }
    }
}

// TODO client-side validation
function run() {
    resultLink.click()
    const formData = new FormData()
    formData.append('fastaText', document.getElementById('fastaText').value);
    selectedGenome = targetGenomes.querySelector('option:checked').value
    formData.append('genome', selectedGenome)
    formData.append('setAmpSize', document.getElementById('setAmpSize').value);
    formData.append('setTmCutoff', document.getElementById('setTmCutoff').value);
    formData.append('setKmer', document.getElementById('setKmer').value);
    formData.append('setDist', document.getElementById('setDist').value);
    formData.append('setCutoffPen', document.getElementById('setCutoffPen').value);
    formData.append('setPenTmDiff', document.getElementById('setPenTmDiff').value);
    formData.append('setPenTmMismatch', document.getElementById('setPenTmMismatch').value);
    formData.append('setPenLength', document.getElementById('setPenLength').value);
    formData.append('setCtmMv', document.getElementById('setCtmMv').value);
    formData.append('setCtmDv', document.getElementById('setCtmDv').value);
    formData.append('setCtmDNA', document.getElementById('setCtmDNA').value);
    formData.append('setCtmDNTP', document.getElementById('setCtmDNTP').value);

    hideElement(resultError)
    showElement(resultInfo)
    hideElement(resultTabs)
    hideElement(sectionResults)

    axios
        .post(`${API_URL}/upload`, formData)
        .then(res => {
            if (res.status === 200) {
                window.data = res.data
                handleSuccess()
            }
        })
        .catch(err => {
            let errorMessage = err
            if (err.response) {
                errorMessage = err.response.data.errors
                .map(error => error.title)
                .join('; ')
            }
            hideElement(resultInfo)
            showElement(resultError)
            resultError.querySelector('#error-message').textContent = errorMessage
        })
}

function handleSuccess() {
    hideElement(resultInfo)
    hideElement(resultError)
    showElement(resultTabs)
    showElement(sectionResults)
    updateResults()
}

function showElement(element) {
  element.classList.remove('d-none')
}

function hideElement(element) {
  element.classList.add('d-none')
}

function checkForUUID() {  
  var path = window.location.search; // .pathname;
  if (path.match(/UUID=.+/)) {
    resultLink.click()

    hideElement(resultError)
    showElement(resultInfo)
    hideElement(resultTabs)
    hideElement(sectionResults)

    var uuid = path.split("UUID=")[1];

    axios
        .get(`${API_URL}/results/` + uuid)
        .then(res => {
            if (res.status === 200) {
                window.data = res.data
                handleSuccess(window.data)
            }
        })
        .catch(err => {
            let errorMessage = err
            if (err.response) {
                errorMessage = err.response.data.errors
                .map(error => error.title)
                .join('; ')
            }
            hideElement(resultInfo)
            showElement(resultError)
            resultError.querySelector('#error-message').textContent = errorMessage
        })
  }
}

function updateResults() {
    var rHTML = ""
    var res = window.data
    var ampcount = res.data.amplicons.length
    var primecount = res.data.primers.length
    var genomeDB = ""
    if (selectedGenome in genomeConv) {
        genomeDB = genomeConv[selectedGenome]
    }
    if (ampcount == 0) {
        rHTML += '<div class="alert alert-warning" role="alert"><strong>No Amplicons Found!</strong></div>\n'
	primerLink.click()
    }
    else {
        rHTML += '<div class="alert alert-success" role="alert"><strong> ' + ampcount + ' Amplicons Found!</strong></div>\n'
    }
    if (primecount == 0) {
        rHTML += '<div class="alert alert-warning" role="alert"><strong>No Primer Binding Sites Found!</strong></div>\n'
    }
    else {
        rHTML += '<div class="alert alert-success" role="alert"><strong> ' + primecount + ' Primer Binding Sites Found!</strong></div>\n'
    }
    rHTML += '<p>Link to this result page:<br />\n'
    rHTML += '<a href="' + `${API_LINK}` + "index.html?UUID=" + res.uuid + '">' + `${API_LINK}` + "index.html?UUID=" + res.uuid + '</a></p>\n'

    sectionResults.innerHTML = '<br />' + rHTML
    rHTML = ""
    if (ampcount > 0) {
        rHTML += '<p><button type="submit" class="btn btn-outline-primary" id="btn-save-amp" onClick="saveAmpAsTSV()">'
        rHTML += '<i class="fas fa-save" style="margin-right: 5px;"></i>Download Amplicons as TSV</button></p>'
        sectionAmpliconLink.innerHTML = rHTML
        rHTML = ""
        if (amStart != 0) {
            showElement(btnAmpUp)
        } else {
            hideElement(btnAmpUp)
        }
        var amStop = amStart + step
        var moreAmp = 1
        if (amStop > ampcount) {
            moreAmp = 0
            amStop = ampcount
        }
        for (var i = amStart; i < amStop; i++) {
            var amp = res.data.amplicons[i]
            rHTML += '<h3>Amplicon ' + (parseInt(amp['Id']) + 1) +'</h3>\n<p>'
            if (genomeDB != "") {
                rHTML += '<a target="_blank" href="' + geneBro + 'db=' + genomeDB + '&position=chr' + amp['Chrom'] + '%3A'
                rHTML += amp['ForPos'] + '%2D' + amp['RevPos'] + '">View region in UCSC Genome Browser</a><br />\n'
            }
            rHTML += '<strong>Length:</strong> ' + amp['Length'] +' bp<br />\n'
            rHTML += '<strong>Penalty:</strong> ' + amp['Penalty'].toFixed(4) +'<br />\n'
            rHTML += '<strong>Location:</strong> ' + amp['Chrom'] + ':' + amp['ForPos'] + '-' + amp['RevPos']+'<br />\n'
            rHTML += '<strong>Forward Primer Name:</strong> ' + amp['ForName'] +'<br />\n'
            rHTML += '<strong>Forward Primer Tm:</strong> ' + amp['ForTm'].toFixed(1) +'&deg;C<br />\n'
            rHTML += '<strong>Forward Primer Sequence:</strong> ' + amp['ForSeq'] +'<br />\n'
            rHTML += '<strong>Reverse Primer Name:</strong> ' + amp['RevName'] +'<br />\n'
            rHTML += '<strong>Reverse Primer Tm:</strong> ' + amp['RevTm'].toFixed(1) +'&deg;C<br />\n'
            rHTML += '<strong>Reverse Primer Sequence:</strong> ' + amp['RevSeq'] +'<br />\n'
            rHTML += '<strong>Amplicon Sequence:</strong><br /><pre>'
            var splitSeq = amp['Seq']
            for (var k = 0; k < splitSeq.length ; k++) {
                if((k % 60 == 0) && (k != 0)) {
                    rHTML += '<br />'
                }
                rHTML += splitSeq.charAt(k)
            }
            rHTML += '</pre><br />\n'
            rHTML += '</p>'
        }
        if (moreAmp == 1) {
            showElement(btnAmpDown)
        } else {
            hideElement(btnAmpDown)
        }
    } else {
        sectionAmpliconLink.innerHTML = ""
    }
    sectionAmpliconData.innerHTML = rHTML
    rHTML = ""
    if (primecount > 0) {
        rHTML += '<p><button type="submit" class="btn btn-outline-primary" id="btn-save-prim" onClick="savePrimAsTSV()">'
        rHTML += '<i class="fas fa-save" style="margin-right: 5px;"></i>Download Primers as TSV</button></p>'
        sectionPrimerLink.innerHTML = rHTML
        rHTML = ""
        if (prStart != 0) {
            showElement(btnPrimerUp)
        } else {
            hideElement(btnPrimerUp)
        }
        var prStop = prStart + step
        var morePrim = 1
        if (prStop > primecount) {
            morePrim = 0
            prStop = primecount
        }
        for (var i = prStart; i < prStop; i++) {
            var prim = res.data.primers[i]
            rHTML += '<h3>Primer Binding Site ' + (parseInt(prim['Id']) + 1) +'</h3>\n<p>'
            if (genomeDB != "") {
                rHTML += '<a target="_blank"href="' + geneBro + 'db=' + genomeDB + '&position=chr' + prim['Chrom'] + '%3A'
                rHTML += prim['Pos'] + '%2D' + (parseInt(prim['Pos']) + prim['Seq'].length) + '">View region in UCSC Genome Browser</a><br />\n'
            }
            rHTML += '<strong>Primer Tm:</strong> ' + prim['Tm'].toFixed(1) +'&deg;C<br />\n'
            rHTML += '<strong>Location:</strong> ' + prim['Chrom'] + ':' + prim['Pos'] + '-' + (parseInt(prim['Pos']) + prim['Seq'].length) + ' on forward<br />\n'
            rHTML += '<strong>Primer Name:</strong> ' + prim['Name'] +'<br />\n'
            rHTML += '<strong>Primer Sequence:</strong> ' + prim['Seq'] +'<br />\n'
            rHTML += '<strong>Genome Sequence:</strong> ' + prim['Genome'] +'<br />\n'
            rHTML += '</p>'
        }
        if (morePrim == 1) {
            showElement(btnPrimerDown)
        } else {
            hideElement(btnPrimerDown)
        }
    } else {
        sectionPrimerLink.innerHTML = ""
    }
    sectionPrimerData.innerHTML = '<br />' + rHTML + '<br />'
}

function ampUp() {
    amStart = amStart - step
    updateResults()
}

function ampDown() {
    amStart = amStart + step
    updateResults()
}

function primerUp() {
    prStart = prStart - step
    updateResults()
}

function primerDown() {
    prStart = prStart + step
    updateResults()
}

window.detectBrowser = detectBrowser;
function detectBrowser() {
    var browser = window.navigator.userAgent.toLowerCase();
    if (browser.indexOf("edge") != -1) {
        return "edge";
    }
    if (browser.indexOf("firefox") != -1) {
        return "firefox";
    }
    if (browser.indexOf("chrome") != -1) {
        return "chrome";
    }
    if (browser.indexOf("safari") != -1) {
        return "safari";
    }
    alert("Unknown Browser: Functionality may be impaired!\n\n" + browser);
    return browser;
}

window.saveJsonFile = saveJsonFile;
function saveJsonFile() {
    if (window.data == "") {
        return;
    }
    var content = JSON.stringify(window.data);
    var a = document.createElement("a");
    document.body.appendChild(a);
    a.style.display = "none";
    var blob = new Blob([content], {type: "application/json"});
    var browser = detectBrowser();
    if (browser != "edge") {
	    var url = window.URL.createObjectURL(blob);
	    a.href = url;
	    a.download = "silica.json";
	    a.click();
	    window.URL.revokeObjectURL(url);
    } else {
        window.navigator.msSaveBlob(blob, fileName);
    }
    return;
};


window.saveAmpAsTSV = saveAmpAsTSV;
function saveAmpAsTSV() {
    if (window.data == "") {
        return;
    }
    var content = "Number\tLength\tPenalty\tLocation\t";
    content += "Forward Name\tForward Tm\tForward Sequence\t";
    content += "Reverse Name\tReverse Tm\tReverse Sequence\t";
    content += "Amplicon Sequence\n";
    var res = window.data
    for (var i = 0; i < res.data.amplicons.length; i++) {
        var amp = res.data.amplicons[i]
        content += (parseInt(amp['Id']) + 1) + "\t" + amp['Length'] + "\t";
        content += amp['Penalty'].toFixed(4) + "\t" + amp['Chrom'] + ':' + amp['ForPos'] + '-' + amp['RevPos'] + "\t";
        content += amp['ForName'] + "\t" + amp['ForTm'].toFixed(1) + "\t" + amp['ForSeq'] + "\t";
        content += amp['RevName'] + "\t" + amp['RevTm'].toFixed(1) + "\t" + amp['RevSeq'] + "\t";
        content += amp['Seq'] + "\n";
    }
    var a = document.createElement("a");
    document.body.appendChild(a);
    a.style.display = "none";
    var blob = new Blob([content], {type: "text/tab-separated-values"});
    var browser = detectBrowser();
    if (browser != "edge") {
	    var url = window.URL.createObjectURL(blob);
	    a.href = url;
	    a.download = "silica_amplicons.tsv";
	    a.click();
	    window.URL.revokeObjectURL(url);
    } else {
        window.navigator.msSaveBlob(blob, fileName);
    }
    return;
};


window.savePrimAsTSV = savePrimAsTSV;
function savePrimAsTSV() {
    if (window.data == "") {
        return;
    }
    var content = "Number\tPrimer Tm\tLocation\tStrand\tPrimer Name\tPrimer Sequence\tGenome Sequence\n";
    var res = window.data
    for (var i = 0; i < res.data.primers.length; i++) {
        var prim = res.data.primers[i]
        content += (parseInt(prim['Id']) + 1) + "\t" + prim['Tm'].toFixed(1) + "\t";
        if (prim['Ori'] == 'reverse') {
            content += prim['Chrom'] + ':' + (parseInt(prim['Pos']) - prim['Seq'].length) + '-' + prim['Pos'] + '\treverse\t'
        } else {
            content += prim['Chrom'] + ':' + prim['Pos'] + '-' + (parseInt(prim['Pos']) + prim['Seq'].length) + '\tforward\t'
        }
        content += prim['Name'] + "\t" + prim['Seq'] + "\t";
        content += prim['Genome'] + "\n";
    }
    var a = document.createElement("a");
    document.body.appendChild(a);
    a.style.display = "none";
    var blob = new Blob([content], {type: "text/tab-separated-values"});
    var browser = detectBrowser();
    if (browser != "edge") {
	    var url = window.URL.createObjectURL(blob);
	    a.href = url;
	    a.download = "silica_primers.tsv";
	    a.click();
	    window.URL.revokeObjectURL(url);
    } else {
        window.navigator.msSaveBlob(blob, fileName);
    }
    return;
};

window.loadJsonFile = loadJsonFile;
function loadJsonFile(f){
    var file = f.target.files[0];
    if (file) { // && file.type.match("text/*")) {
        var reader = new FileReader();
        reader.onload = function(event) {
            window.data = JSON.parse(event.target.result);
            handleSuccess();
        }
        reader.readAsText(file);
    } else {
        alert("Error opening file");
    }
}
