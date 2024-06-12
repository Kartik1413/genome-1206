function showInput(type) {
  // Hide all input fields and buttons
  document.getElementById('input-smiles').style.display = 'none';
  document.getElementById('btn-predict-smiles').style.display = 'none';
  document.getElementById('input-molecular').style.display = 'none';
  document.getElementById('btn-predict-molecular').style.display = 'none';

  // Clear result tables
  document.getElementById('table-body').innerHTML = '';
  document.getElementById('result-table').style.display = 'none';
  document.getElementById('table-body-molecular').innerHTML = '';
  document.getElementById('result-table-molecular').style.display = 'none';

  // Show the selected input field and button
  if (type === 'smiles') {
     document.getElementById('input-smiles').style.display = 'block';
     document.getElementById('btn-predict-smiles').style.display = 'block';
  } else if (type === 'molecular') {
     document.getElementById('input-molecular').style.display = 'block';
     document.getElementById('btn-predict-molecular').style.display = 'block';
  }
}

function openModal(modalId) {
  var modal = document.getElementById(modalId);
  modal.style.display = "block";
  resetModalInputs(modalId);
}

function closeModal(modalId) {
  var modal = document.getElementById(modalId);
  modal.style.display = "none";
  resetModalInputs(modalId);
}

function resetModalInputs(modalId) {
  if (modalId === 'modal_single_smile') {
     document.getElementById('input-smiles').value = '';
     document.getElementById('table-body').innerHTML = '';
     document.getElementById('result-table').style.display = 'none';
     document.getElementById('input-molecular').value = '';
     document.getElementById('table-body-molecular').innerHTML = '';
     document.getElementById('result-table-molecular').style.display = 'none';
  }
}

function predictSmiles() {
  const smiles = document.getElementById('input-smiles').value;
  fetch('/predict1', {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json'
    },
    body: JSON.stringify({ smiles: smiles })
  })
  .then(response => response.json())
  .then(data => {
    const tableBody = document.getElementById('table-body');
    tableBody.innerHTML = '';

    for (const [key, value] of Object.entries(data.result)) {
      const row = document.createElement('tr');
      const keyCell = document.createElement('td');
      const valueCell = document.createElement('td');
      keyCell.textContent = key;
      valueCell.textContent = value;
      row.appendChild(keyCell);
      row.appendChild(valueCell);
      tableBody.appendChild(row);
    }

    document.getElementById('result-table').style.display = 'table';
  })
  .catch((error) => {
    console.error('Error:', error);
  });
}

function predictMolecular() {
  const formula = document.getElementById('input-molecular').value;
  fetch('/predict2', {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json'
    },
    body: JSON.stringify({ formula: formula })
  })
  .then(response => response.json())
  .then(data => {
    const tableBody = document.getElementById('table-body-molecular');
    tableBody.innerHTML = '';

    for (const [key, value] of Object.entries(data.result)) {
      const row = document.createElement('tr');
      const keyCell = document.createElement('td');
      const valueCell = document.createElement('td');
      keyCell.textContent = key;
      valueCell.textContent = value;
      row.appendChild(keyCell);
      row.appendChild(valueCell);
      tableBody.appendChild(row);
    }

    document.getElementById('result-table-molecular').style.display = 'table';
  })
  .catch((error) => {
    console.error('Error:', error);
  });
}
