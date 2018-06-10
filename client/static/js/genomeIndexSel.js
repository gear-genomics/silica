/* global XMLHttpRequest */

var genBtn = document.getElementById('genomeIndexSel')
genBtn.addEventListener('load', buildGenomeIndexButton)

function buildGenomeIndexButton() {
    var loca = 'http://0.0.0.0:3300';
    if (location.origin.startsWith("http")) {
        loca = location.origin;
    }
    var req = new XMLHttpRequest()
    req.addEventListener('load', bgiResults)
    req.open('GET', loca + '/genomeindex', true)
    req.send()
}

function bgiResults() {
    if (this.status === 200) {
        bgiProcess(this.response)
    } else {
        alert("Fatal error loading genome index!") 
    }
}

function bgiProcess(data) {
    var res = JSON.parse(data)
    var rhtml = '<input class="span2" id="genome" name="genome" type="hidden">\n<div class="alert alert-primary mb-2" role="alert">\n'
    rhtml += 'and select genome\n<br />\n<p></p>\n<div class="dropdown">\n'
    rhtml += '<button type="button" class="btn btn-secondary dropdown-toggle" data-toggle="dropdown" aria-haspopup="true"'
    rhtml += ' aria-expanded="false" id="btn-genome" name="btn-genome">\nGenome\n</button>\n<div class="dropdown-menu">'
    for (var i = 0; i < res.length; i++) {
        rhtml += '<a class="dropdown-item" href="#" onclick="bgiUp(\'' + res[i].name + '\',\'' + res[i].file + '\')">' + res[i].name + '</a>'
    }
    rhtml += '</div>\n</div>\n</div>\n</div>\n'
    genBtn.innerHTML = rhtml
}

function bgiUp(gname,gfile) {
    document.getElementById('genome').value = gfile;
    document.getElementById('btn-genome').innerHTML = gname;
}



