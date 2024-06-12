document.addEventListener('DOMContentLoaded', function () {
    var generateButton = document.getElementById('generate-button');
    if (generateButton) {
        generateButton.addEventListener('click', function () {
            uploadFile();
        });
    }
});

function uploadFile() {
    var fileInput = document.getElementById('file-input-pubchem');
    var files = fileInput.files;

    if (files.length === 0) {
        alert('Please select a file to upload.');
        return;
    }

    var formData = new FormData();
    formData.append('file', files[0]);

    var xhr = new XMLHttpRequest();
    xhr.open('POST', '/upload');
    xhr.onload = function () {
        if (xhr.status === 200) {
            var response = JSON.parse(xhr.responseText);
            displayContents(response.contents);
        } else {
            alert('Error uploading file. Please try again.');
        }
    };
    xhr.send(formData);
}

function displayContents(contents) {
    var table = document.getElementById('file-contents');
    table.innerHTML = ''; // Clear previous contents

    // Create table header row
    var thead = document.createElement('thead');
    var headerRow = document.createElement('tr');
    Object.keys(contents[0]).forEach(function (key) {
        var th = document.createElement('th');
        th.textContent = key;
        headerRow.appendChild(th);
    });
    thead.appendChild(headerRow);
    table.appendChild(thead);

    // Create table body rows
    var tbody = document.createElement('tbody');
    contents.forEach(function (row) {
        var tr = document.createElement('tr');
        Object.values(row).forEach(function (cellContent) {
            var td = document.createElement('td');
            td.textContent = cellContent;
            tr.appendChild(td);
        });
        tbody.appendChild(tr);
    });
    table.appendChild(tbody);
}

function closeModal(modalId) {
    var modal = document.getElementById(modalId);
    modal.style.display = "none";
}
