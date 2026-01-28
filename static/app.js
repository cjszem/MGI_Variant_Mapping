// Triggered when user clicks "Submit"
document.getElementById("runButton").addEventListener("click", function () {
    const text = document.getElementById("variantInput").value;

    fetch("http://localhost:8000/run_variants", {
        method: "POST",
        headers: {"Content-Type": "application/json"},
        body: JSON.stringify({ input: text })
    })
    .then(res => res.json())
    .then(data => {
        // store full result globally so we can: (1) filter per input, (2) download raw JSON
        window.fullResults = data;

        // optional: show raw JSON somewhere while debugging
        document.getElementById("output").textContent = JSON.stringify(data, null, 2);

        buildSidebar(data);

        // auto-select first entry if present
        if (data.human_genes && data.human_genes.length > 0) {
            renderTablesForIndex(0, data);
        }
    })
    .catch(err => {
        console.error(err);
        document.getElementById("output").textContent = "Error running query";
    });
});


// Build the sidebar with one button per input variant
function buildSidebar(data) {
    const sidebar = document.getElementById("sidebar");
    sidebar.innerHTML = "<h3>Inputs</h3>";

    if (!data.human_genes || data.human_genes.length === 0) {
        sidebar.innerHTML += "<p>No inputs found.</p>";
        return;
    }

    data.human_genes.forEach((row, idx) => {
        const btn = document.createElement("button");
        btn.className = "sidebar-button";
        btn.textContent = `${row.Gene} (${row.Chr}:${row.Start}-${row.End})`;
        btn.onclick = () => renderTablesForIndex(idx, data);
        sidebar.appendChild(btn);
    });
}


// Helper: turn an array of objects into an HTML table
function jsonToHTMLTable(json) {
    if (!json || json.length === 0) {
        return "<p>No data available.</p>";
    }

    const columns = Object.keys(json[0]);
    const header = columns.map(c => `<th>${c}</th>`).join("");

    const rows = json.map(row =>
        `<tr>${columns.map(c => `<td>${row[c] ?? ""}</td>`).join("")}</tr>`
    ).join("");

    return `
        <table class="result-table">
            <thead><tr>${header}</tr></thead>
            <tbody>${rows}</tbody>
        </table>
    `;
}


// Helper: filter by gene + start + end for a given input
function filterByIndex(df, gene, start, end) {
    if (!df) return [];
    return df.filter(row =>
        row.Gene === gene &&
        Number(row.Start) === Number(start) &&
        Number(row.End) === Number(end)
    );
}


// Render all tables for the selected input
function renderTablesForIndex(idx, data) {
    const entry = data.human_genes[idx];

    const humGenes = filterByIndex(data.human_genes, entry.Gene, entry.Start, entry.End);
    const humProteins = filterByIndex(data.human_proteins, entry.Gene, entry.Start, entry.End);
    const mouseGenes = filterByIndex(data.mouse_genes, entry.Gene, entry.Start, entry.End);
    const mouseProteins = filterByIndex(data.mouse_proteins, entry.Gene, entry.Start, entry.End);
    const scores = data.scores || [];

    const output = document.getElementById("output");
    output.innerHTML = `
        <h2>Human Gene</h2>
        ${jsonToHTMLTable(humGenes)}

        <h2>Human Proteins</h2>
        ${jsonToHTMLTable(humProteins)}

        <h2>Mouse Gene</h2>
        ${jsonToHTMLTable(mouseGenes)}

        <h2>Mouse Proteins</h2>
        ${jsonToHTMLTable(mouseProteins)}

        <h2>Similarity Scores</h2>
        ${jsonToHTMLTable(scores)}
    `;
}


// Download full JSON for all data
document.getElementById("downloadJSON").addEventListener("click", function () {
    if (!window.fullResults) {
        alert("No results to download yet.");
        return;
    }

    const dataStr = "data:text/json;charset=utf-8," +
        encodeURIComponent(JSON.stringify(window.fullResults, null, 2));

    const a = document.createElement("a");
    a.setAttribute("href", dataStr);
    a.setAttribute("download", "variant_mapping_results.json");
    a.click();
});