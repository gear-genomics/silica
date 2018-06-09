/* global XMLHttpRequest */

var res = ""
var prStart = 0
var amStart = 0
var step = 30

var fileLoad = document.getElementById('fasta');
fileLoad.addEventListener('change', loadFasta, false);
var submitButton = document.getElementById('btn-submit')
submitButton.addEventListener('click', submit)
var sampleButton = document.getElementById('btn-example')
sampleButton.addEventListener('click', sampleData)
var helpButton = document.getElementById('btn-help')
helpButton.addEventListener('click', goToHelp)
const resultLink = document.getElementById('link-results')
const helpLink = document.getElementById('link-help')
const primerLink = document.getElementById('link-primers')
var sectionResults = document.getElementById('results')
var sectionAmplicons = document.getElementById('res-amplicons')
var sectionPrimers = document.getElementById('res-primers')

$('#mainTab a').on('click', function (e) {
  e.preventDefault()
  $(this).tab('show')
})

$('#resTab a').on('click', function (e) {
  e.preventDefault()
  $(this).tab('show')
})

var spinnerHtml = '<p>Analysis takes a couple of seconds to run, please be patient.</p><i class="fa fa-spinner fa-pulse fa-3x fa-fw"></i><br /><br />'

function checkPath(){
    var path = window.location.pathname
    if (path != "/") {
        loadLink (path)
    }
}

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
    document.getElementById('genome').value = 'Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz';
    document.getElementById('btn-genome').innerHTML = 'Homo Sapiens GRCh37';
}

function submit () {
    var data = new FormData();
    data.append('fastaText', document.getElementById('fastaText').value);
    data.append('genome', document.getElementById('genome').value);
    data.append('setAmpSize', document.getElementById('setAmpSize').value);
    data.append('setTmCutoff', document.getElementById('setTmCutoff').value);
    data.append('setKmer', document.getElementById('setKmer').value);
    data.append('setDist', document.getElementById('setDist').value);
    data.append('setCutoffPen', document.getElementById('setCutoffPen').value);
    data.append('setPenTmDiff', document.getElementById('setPenTmDiff').value);
    data.append('setPenTmMismatch', document.getElementById('setPenTmMismatch').value);
    data.append('setPenLength', document.getElementById('setPenLength').value);
    data.append('setCtmMv', document.getElementById('setCtmMv').value);
    data.append('setCtmDv', document.getElementById('setCtmDv').value);
    data.append('setCtmDNA', document.getElementById('setCtmDNA').value);
    data.append('setCtmDNTP', document.getElementById('setCtmDNTP').value);
    doSubmit (data);
}

function loadLink (uuid) {
    resultLink.click();
    var loca = 'http://0.0.0.0:3300';
    if (location.origin.startsWith("http")) {
        loca = location.origin;
    }
    var req = new XMLHttpRequest()
    req.addEventListener('load', displayResults)
    req.open('GET', loca + '/results' + uuid, true)
    req.send()
    sectionResults.innerHTML = spinnerHtml
}

function doSubmit (data) {
    resultLink.click();
    var loca = 'http://0.0.0.0:3300';
    if (location.origin.startsWith("http")) {
        loca = location.origin;
    }
    var req = new XMLHttpRequest()
    req.addEventListener('load', displayResults)
    req.open('POST', loca + '/upload', true)
    req.send(data)
    sectionResults.innerHTML = spinnerHtml
}

function displayResults() {
    if (this.status === 200) {
        displayData(this.response)
    } else {
        displayError(this.response)
    }
    resultLink.click();
}

function displayData(data) {
    res = JSON.parse(data)
    updateResults()
}

function displayError(data) {
    alert(data)
    var res = JSON.parse(data)
    for (var i = 0; i < res["errors"].length; i++) {
        sectionResults.innerHTML = '<br /><div class="error">' + res["errors"][i]['title'] + '</div><br />'
    }
    sectionAmplicons.innerHTML = '<br />'
    sectionPrimers.innerHTML = '<br />'
}

function updateResults() {
    var rHTML = ""
    var ampcount = res.data.amplicon.length
    var primecount = res.data.primer.length
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
    sectionResults.innerHTML = '<br />' + rHTML + '<br />'
    rHTML = ""
    if (ampcount > 0) {
        rHTML += '<p>Download all as <a href="/silica/download' + res.data.link + '-ac">CSV</a>'
        rHTML += ' or <a href="/silica/download/' + res.data.link + '-aj">JSON</a></p>'
        if (amStart != 0) {
            rHTML += '<button type="button" class="btn btn-secondary" onclick="ampUp()">Up</button>'
        }
        var amStop = amStart + step
        var moreAmp = 1
        if (amStop > ampcount) {
            moreAmp = 0
            amStop = ampcount
        }
        for (var i = amStart; i < amStop; i++) {
            var amp = res.data.amplicon[i]
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
            rHTML += '<button type="button" class="btn btn-secondary" onclick="ampDown()">Down</button>'
        }
    }
    sectionAmplicons.innerHTML = '<br />' + rHTML + '<br />'
    rHTML = ""
    if (primecount > 0) {
        rHTML += '<p>Download all as <a href="/silica/download' + res.data.link + '-pc">CSV</a>'
        rHTML += ' or <a href="/silica/download/' + res.data.link + '-pj">JSON</a></p>'
        if (prStart != 0) {
            rHTML += '<button type="button" class="btn btn-secondary" onclick="primerUp()">Up</button>'
        }
        var prStop = prStart + step
        var morePrim = 1
        if (prStop > primecount) {
            morePrim = 0
            prStop = primecount
        }
        for (var i = prStart; i < prStop; i++) {
            var prim = res.data.primer[i]
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
            rHTML += '<button type="button" class="btn btn-secondary" onclick="primerDown()">Down</button>'
	}
    }
    sectionPrimers.innerHTML = '<br />' + rHTML + '<br />'
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

