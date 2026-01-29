// Run variant-model mapping for submitted variants on click.
document.getElementById('runButton').addEventListener('click', function () {
    const text = document.getElementById('variantInput').value;

    fetch('http://localhost:8000/run_variants', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({ input: text })
    })
    .then(res => res.json())
    .then(data => {
        // store full results
        window.fullJSON = data;

        // build sidebar
        buildSidebar(data);
    })
    .catch(err => {console.error(err);});
});


// Build the sidebar with one button per input variant
function buildSidebar(data) {
    const sidebar = document.getElementById('sidebar');

    // Create sidebar header
    sidebar.innerHTML = '<h3>Inputs:</h3>';

    // Handle no results
    if (!data.human_genes || data.human_genes.length === 0) {
        sidebar.innerHTML += '<h2>No inputs found.</h2>';
        return;
    }

    // Create sidebar buttons
    data.human_proteins.forEach((row, idx) => {
        const btn = document.createElement('button');
        btn.className = 'sidebar-button';
        btn.textContent = row['Input'];
        btn.onclick = () => {
            // deactivate all other buttons
            document.querySelectorAll('#sidebar .sidebar-button')
            .forEach(b => b.classList.remove('active'));

            // activate clicked button
            btn.classList.add('active')

            // render output tables
            renderTables(row['Input'], row['Gene Symbol']); 
        }
        sidebar.appendChild(btn);
    });

    // Unhide results
    document.getElementById("layoutOutput").style.display = "flex";
}


// Render all tables for the selected input
function renderTables(input, gene) {

    // Extract homolog
    homolog = window.fullJSON.gene_mapping.find(r => r['Hum Gene'] === gene)?.['Mus Gene']

    // Extract relevant tables
    const humGenes = window.fullJSON.human_genes.filter(row => row['Gene Symbol'] === gene);
    const humProteins = window.fullJSON.human_proteins.filter(row => row['Input'] === input);
    const mouseGenes = window.fullJSON.mouse_genes.filter(row => row['Gene Symbol'] === homolog);
    const mouseProteins = window.fullJSON.mouse_proteins.filter(row => row['Gene Symbol'] === homolog);
    const scores = window.fullJSON.scores.filter(row => row['Input'] === input);
    const phenotypes = window.fullJSON.phenotypes.filter(row => row['Gene Symbol'] === homolog);

    // Create HTML tables
    results.innerHTML = `
        <div class='tableWrapper'>
            <h2>Human Gene</h2>
            ${jsonToHTMLTable(humGenes, "humGenesTable")}
        </div>

        <div class='tableWrapper'>
            <h2>Human Variant</h2>
            ${jsonToHTMLTable(humProteins, "humVariantTable")}
        </div>

        <div class='tableWrapper'>
            <h2>Mouse Homolog</h2>
            ${jsonToHTMLTable(mouseGenes, "musGenesTable")}
        </div>

        <div class='tableWrapper'>
            <h2>Mouse Models</h2>
            ${jsonToHTMLTable(mouseProteins, "musModelTable")}
        </div>

        <div class='tableWrapper'>
            <h2>Similarity Scores</h2>
            ${jsonToHTMLTable(scores, "scoresTable")}
        </div>


        <div class='tableWrapper'>
            <h2>Mouse Allele Phenotypes</h2>
            ${jsonToHTMLTable(phenotypes, "phenotypeTable")}
        </div>
    `;


    // Create DataTable elements
    $('#humGenesTable').DataTable({
        columnDefs: [{ targets: [], visible: false }],
        dom: 't',
        pageLength: 20,
    });

    $('#humVariantTable').DataTable({
        columnDefs: [{ targets: [0, 1, 13, 14, 17], visible: false }], // Hide input, submission, refAA, varAA, and MondoSet
        dom: 't',
        pageLength: 20,
    });

    $('#musGenesTable').DataTable({
        columnDefs: [{ targets: [], visible: false }],
        dom: 't',
        pageLength: 20,
    });

    $('#musModelTable').DataTable({
        columnDefs: [{ targets: [11, 12, 15], visible: false }], // Hide refAA, varAA, and MondoSet
        dom: 'tp',
        pageLength: 20,
    });

    $('#scoresTable').DataTable({
        columnDefs: [{ targets: [0], visible: false }], // Hide input
        dom: 'tp',
        pageLength: 20,
    });

    $('#phenotypeTable').DataTable({
        columnDefs: [{ targets: [], visible: false }],
        dom: 'tp',
        pageLength: 20,
    });
}


// Turn an array of objects into an HTML table
function jsonToHTMLTable(data, tableId) {
    if (!data || data.length === 0) return "<p>No data</p>";

    const columns = Object.keys(data[0]);

    let html = `<table id="${tableId}" class="datatable">`;
    html += "<thead><tr>";

    columns.forEach(col => {
        html += `<th>${col}</th>`;
    });

    html += "</tr></thead><tbody>";

    data.forEach(row => {
        html += "<tr>";
        columns.forEach(col => {
            html += `<td>${row[col] !== undefined ? row[col] : ""}</td>`;
        });
        html += "</tr>";
    });

    html += "</tbody></table>";

    return html;
}


// Download full JSON for all data
document.getElementById('downloadJSON').addEventListener('click', function () {
    if (!window.fullJSON) {
        alert('No results yet.');
        return;
    }

    const dataStr = 'data:text/json;charset=utf-8,' +
        encodeURIComponent(JSON.stringify(window.fullJSON, null, 2));

    const a = document.createElement('a');
    a.setAttribute('href', dataStr);
    a.setAttribute('download', 'variant_mapping_results.json');
    a.click();
});