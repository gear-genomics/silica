const API_URL = process.env.API_URL
const API_LINK = process.env.API_LINK

var res = ""
var prStart = 0
var amStart = 0
var step = 30

var fileLoad = document.getElementById('fasta');
fileLoad.addEventListener('change', loadFasta, false);
var submitButton = document.getElementById('btn-submit')
submitButton.addEventListener('click', run)
var sampleButton = document.getElementById('btn-example')
sampleButton.addEventListener('click', sampleData)
var helpButton = document.getElementById('btn-help')
helpButton.addEventListener('click', goToHelp)


const resultLink = document.getElementById('link-results')
const helpLink = document.getElementById('link-help')
const resultInfo = document.getElementById('result-info')
const resultError = document.getElementById('result-error')
const resultTabs = document.getElementById('result-tabs')
const primerLink = document.getElementById('link-primers')
const targetGenomes = document.getElementById('target-genome')
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
    const genome = targetGenomes.querySelector('option:checked').value
    formData.append('genome', genome)
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
                handleSuccess(res.data)
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

function handleSuccess(data) {
    hideElement(resultInfo)
    hideElement(resultError)
    showElement(resultTabs)
    showElement(sectionResults)
    res = data
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
                handleSuccess(res.data)
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
    var ampcount = res.data.amplicons.length
    var primecount = res.data.primers.length
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
    rHTML += '<a href="' + `${API_LINK}` + "index.html?UUID=" + res.data.uuid + '">' + `${API_LINK}` + "index.html?UUID=" + res.data.uuid + '</a></p>\n'
    sectionResults.innerHTML = '<br />' + rHTML + '<br />'
    rHTML = ""
    if (ampcount > 0) {
        rHTML += '<p>Download all as <a href="' + `${API_URL}/download/` + res.data.uuid + '-ac" download>CSV</a>'
        rHTML += ' or <a href="' + `${API_URL}/download/` + res.data.uuid + '-aj" download>JSON</a></p>'
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
            rHTML += '<strong>Length:</strong> ' + amp['Length'] +' bp<br />\n'
            rHTML += '<strong>Penalty:</strong> ' + amp['Penalty'] +'<br />\n'
            rHTML += '<strong>Location:</strong> ' + amp['Chrom'] + ':' + amp['ForPos'] + '-' + amp['RevPos']+'<br />\n'
            rHTML += '<strong>Forward Primer Name:</strong> ' + amp['ForName'] +'<br />\n'
            rHTML += '<strong>Forward Primer Tm:</strong> ' + amp['ForTm'] +'&deg;C<br />\n'
            rHTML += '<strong>Forward Primer Sequence:</strong> ' + amp['ForSeq'] +'<br />\n'
            rHTML += '<strong>Reverse Primer Name:</strong> ' + amp['RevName'] +'<br />\n'
            rHTML += '<strong>Reverse Primer Tm:</strong> ' + amp['RevTm'] +'&deg;C<br />\n'
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
        rHTML += '<p>Download all as <a href="' + `${API_URL}/download/` + res.data.uuid + '-pc" download>CSV</a>'
        rHTML += ' or <a href="' + `${API_URL}/download/` + res.data.uuid + '-pj" download>JSON</a></p>'
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
            rHTML += '<strong>Primer Tm:</strong> ' + prim['Tm'] +'&deg;C<br />\n'
            if (prim['Ori'] == 'reverse') {
                rHTML += '<strong>Location:</strong> ' + prim['Chrom'] + ':' + (parseInt(prim['Pos']) - prim['Seq'].length) + '-' + prim['Pos'] + ' on reverse<br />\n'
            } else {
                rHTML += '<strong>Location:</strong> ' + prim['Chrom'] + ':' + prim['Pos'] + '-' + (parseInt(prim['Pos']) + prim['Seq'].length) + ' on forward<br />\n'
            }
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

