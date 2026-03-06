// API Configuration
const API_BASE_URL = '/api/v1';

// Example variants (hg38/GRCh38 coordinates - automatic conversion to hg19 when needed)
const examples = {
    missense: {
        gene: 'APOE',
        chromosome: '19',
        position: '45411941',
        reference: 'T',
        alternate: 'C',
        hgvs_c: 'c.388T>C',
        hgvs_p: 'p.Cys130Arg'
    },
    frameshift: {
        gene: 'BRCA1',
        chromosome: '17',
        position: '43124027',
        reference: 'C',
        alternate: 'CA',
        hgvs_c: 'c.5266dup',
        hgvs_p: 'p.Gln1756fs'
    },
    splicing: {
        gene: 'CFTR',
        chromosome: '7',
        position: '117559590',
        reference: 'ATC',
        alternate: 'A',
        hgvs_c: 'c.1521_1523delCTT',
        hgvs_p: 'p.Phe508del'
    },
    synonymous: {
        gene: 'BRCA1',
        chromosome: '17',
        position: '43091825',
        reference: 'A',
        alternate: 'G',
        hgvs_c: 'c.4308T>C',
        hgvs_p: 'p.(=)'
    },
    utr5: {
        gene: 'CDKN2A',
        chromosome: '9',
        position: '21974826',
        reference: 'G',
        alternate: 'A',
        hgvs_c: 'c.-34G>A',
        hgvs_p: ''
    },
    utr3: {
        gene: 'PTEN',
        chromosome: '10',
        position: '89692910',
        reference: 'C',
        alternate: 'T',
        hgvs_c: 'c.*95C>T',
        hgvs_p: ''
    }
};

// Load example variant
function loadExample(type) {
    const example = examples[type];
    if (!example) return;

    document.getElementById('gene').value = example.gene;
    document.getElementById('chromosome').value = example.chromosome;
    document.getElementById('position').value = example.position;
    document.getElementById('reference').value = example.reference;
    document.getElementById('alternate').value = example.alternate;
    document.getElementById('hgvs_c').value = example.hgvs_c || '';
    document.getElementById('hgvs_p').value = example.hgvs_p || '';
}

// Form submission handler
document.getElementById('variantForm').addEventListener('submit', async (e) => {
    e.preventDefault();

    // Get form data
    const variantData = {
        gene: document.getElementById('gene').value.trim(),
        chromosome: document.getElementById('chromosome').value.trim().replace('chr', ''),
        position: parseInt(document.getElementById('position').value),
        reference: document.getElementById('reference').value.trim().toUpperCase(),
        alternate: document.getElementById('alternate').value.trim().toUpperCase(),
        hgvs_c: document.getElementById('hgvs_c').value.trim() || null,
        hgvs_p: document.getElementById('hgvs_p').value.trim() || null
    };

    // Show loading spinner
    showLoading();

    try {
        // Call API
        const response = await fetch(`${API_BASE_URL}/analyze`, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify(variantData)
        });

        if (!response.ok) {
            throw new Error(`API Error: ${response.status} ${response.statusText}`);
        }

        const result = await response.json();

        // Display results
        displayResults(result);

    } catch (error) {
        console.error('Error:', error);
        showError(`Errore nell'analisi: ${error.message}`);
    } finally {
        hideLoading();
    }
});

function showLoading() {
    document.getElementById('loadingSpinner').style.display = 'block';
    document.getElementById('resultsContainer').style.display = 'none';
    document.getElementById('errorMessage').style.display = 'none';
}

function hideLoading() {
    document.getElementById('loadingSpinner').style.display = 'none';
}

function showError(message) {
    const errorDiv = document.getElementById('errorMessage');
    errorDiv.textContent = message;
    errorDiv.style.display = 'block';
    document.getElementById('resultsContainer').style.display = 'none';
}

function displayResults(data) {
    // Show results container
    document.getElementById('resultsContainer').style.display = 'block';
    document.getElementById('errorMessage').style.display = 'none';

    // Display classification summary
    displayClassificationSummary(data);

    // Display variant information with population data
    displayVariantInfo(data.variant, data.population_data);

    // Display predictions
    displayPredictions(data.predictions, data.population_data);

    // Display literature
    displayLiterature(data.literature);

    // Display ACMG classification
    displayACMG(data.acmg_classification);

    // Display external links
    displayExternalLinks(data.tool_links);

    // Scroll to results
    document.getElementById('resultsContainer').scrollIntoView({ behavior: 'smooth' });
}

function displayClassificationSummary(data) {
    const container = document.getElementById('classificationSummary');
    const classification = data.acmg_classification?.classification || 'Unknown';
    const summary = data.summary || {};

    const classificationClass = classification
        .toLowerCase()
        .replace(/\s+/g, '-')
        .replace(/[()]/g, '');

    let html = `
        <div class="classification-badge classification-${classificationClass}">
            ${classification}
        </div>
        <div style="margin-top: 20px;">
            <h3>Riepilogo</h3>
    `;

    if (summary.key_findings && summary.key_findings.length > 0) {
        html += '<h4 style="margin-top: 15px; color: var(--text-secondary);">Risultati Principali:</h4><ul>';
        summary.key_findings.forEach(finding => {
            html += `<li>${finding}</li>`;
        });
        html += '</ul>';
    }

    if (summary.recommendations && summary.recommendations.length > 0) {
        html += '<h4 style="margin-top: 15px; color: var(--text-secondary);">Raccomandazioni:</h4><ul>';
        summary.recommendations.forEach(rec => {
            html += `<li>${rec}</li>`;
        });
        html += '</ul>';
    }

    html += '</div>';
    container.innerHTML = html;
}

function displayVariantInfo(variant, populationData) {
    const container = document.getElementById('variantInfo');
    const popData = populationData || {};

    let html = `
        <div class="info-grid">
            <div class="info-item">
                <div class="info-label">Gene</div>
                <div class="info-value">${variant.gene || 'N/A'}</div>
            </div>
            <div class="info-item">
                <div class="info-label">Posizione Genomica</div>
                <div class="info-value">chr${variant.chromosome}:${variant.position}</div>
            </div>
            <div class="info-item">
                <div class="info-label">Cambiamento</div>
                <div class="info-value">${variant.reference} → ${variant.alternate}</div>
            </div>
            <div class="info-item">
                <div class="info-label">Tipo Variante</div>
                <div class="info-value"><span class="variant-type">${variant.variant_type || 'Unknown'}</span></div>
            </div>
    `;

    // dbSNP rsID
    if (popData.dbsnp_rsid) {
        html += `
            <div class="info-item">
                <div class="info-label">dbSNP ID</div>
                <div class="info-value">
                    <a href="${popData.dbsnp_url}" target="_blank" style="color: var(--primary-color); font-weight: 600;">
                        ${popData.dbsnp_rsid}
                    </a>
                </div>
            </div>
        `;
    }

    // gnomAD Frequencies
    if (popData.gnomad_max_af !== null && popData.gnomad_max_af !== undefined) {
        const af = popData.gnomad_max_af;
        const afFormatted = af < 0.0001 ? af.toExponential(2) : af.toFixed(6);
        html += `
            <div class="info-item">
                <div class="info-label">gnomAD Freq (Max)</div>
                <div class="info-value" style="font-weight: 600; color: ${af > 0.01 ? 'var(--danger-color)' : af > 0.001 ? 'var(--warning-color)' : 'var(--success-color)'};">
                    ${afFormatted}
                </div>
            </div>
        `;
    }

    // ClinVar Significance
    if (popData.clinvar_significance) {
        html += `
            <div class="info-item">
                <div class="info-label">ClinVar</div>
                <div class="info-value">
                    <a href="${popData.clinvar_url}" target="_blank" style="color: var(--primary-color);">
                        ${popData.clinvar_significance}
                    </a>
                </div>
            </div>
        `;
    }

    if (variant.hgvs_c) {
        html += `
            <div class="info-item">
                <div class="info-label">HGVS Coding</div>
                <div class="info-value">${variant.hgvs_c}</div>
            </div>
        `;
    }

    if (variant.hgvs_p) {
        html += `
            <div class="info-item">
                <div class="info-label">HGVS Protein</div>
                <div class="info-value">${variant.hgvs_p}</div>
            </div>
        `;
    }

    html += '</div>';
    container.innerHTML = html;
}

function displayPredictions(predictions, populationData) {
    const container = document.getElementById('predictions');
    const popData = populationData || {};
    let html = '';

    // Real prediction scores from MyVariant (dbNSFP + CADD)
    if (popData.found && (popData.sift_score || popData.polyphen2_score || popData.cadd_phred || popData.revel_score)) {
        html += '<h3>🎯 Predizioni da Database (Valori Reali)</h3><div class="prediction-grid">';

        // SIFT
        if (popData.sift_score !== null && popData.sift_score !== undefined) {
            const siftColor = popData.sift_score <= 0.05 ? 'var(--danger-color)' : 'var(--success-color)';
            html += `
                <div class="prediction-item" style="border-left-color: ${siftColor};">
                    <h4>SIFT</h4>
                    <div class="prediction-score" style="color: ${siftColor};">${popData.sift_score.toFixed(3)}</div>
                    <div><strong>Predizione:</strong> ${popData.sift_prediction || (popData.sift_score <= 0.05 ? 'Deleterious' : 'Tolerated')}</div>
                    <div style="margin-top: 8px; font-size: 0.9rem;">≤0.05 = dannosa</div>
                    <div style="margin-top: 10px;"><a href="https://sift.bii.a-star.edu.sg/" target="_blank" style="color: var(--primary-color);">Vai al tool →</a></div>
                </div>
            `;
        }

        // PolyPhen-2
        if (popData.polyphen2_score !== null && popData.polyphen2_score !== undefined) {
            const polyphenColor = popData.polyphen2_score >= 0.85 ? 'var(--danger-color)' : popData.polyphen2_score >= 0.15 ? 'var(--warning-color)' : 'var(--success-color)';
            html += `
                <div class="prediction-item" style="border-left-color: ${polyphenColor};">
                    <h4>PolyPhen-2</h4>
                    <div class="prediction-score" style="color: ${polyphenColor};">${popData.polyphen2_score.toFixed(3)}</div>
                    <div><strong>Predizione:</strong> ${popData.polyphen2_prediction || (popData.polyphen2_score >= 0.85 ? 'Probably damaging' : popData.polyphen2_score >= 0.15 ? 'Possibly damaging' : 'Benign')}</div>
                    <div style="margin-top: 8px; font-size: 0.9rem;">>0.85 = probably damaging</div>
                    <div style="margin-top: 10px;"><a href="http://genetics.bwh.harvard.edu/pph2/" target="_blank" style="color: var(--primary-color);">Vai al tool →</a></div>
                </div>
            `;
        }

        // CADD
        if (popData.cadd_phred !== null && popData.cadd_phred !== undefined) {
            const caddColor = popData.cadd_phred >= 30 ? 'var(--danger-color)' : popData.cadd_phred >= 20 ? 'var(--warning-color)' : 'var(--success-color)';
            html += `
                <div class="prediction-item" style="border-left-color: ${caddColor};">
                    <h4>CADD</h4>
                    <div class="prediction-score" style="color: ${caddColor};">PHRED: ${popData.cadd_phred.toFixed(1)}</div>
                    ${popData.cadd_raw ? `<div>Raw: ${popData.cadd_raw.toFixed(3)}</div>` : ''}
                    <div style="margin-top: 8px; font-size: 0.9rem;">${popData.cadd_interpretation || '>20 = top 1% deleterious'}</div>
                    <div style="margin-top: 10px;"><a href="https://cadd.gs.washington.edu/" target="_blank" style="color: var(--primary-color);">Vai al tool →</a></div>
                </div>
            `;
        }

        // REVEL
        if (popData.revel_score !== null && popData.revel_score !== undefined) {
            const revelColor = popData.revel_score >= 0.75 ? 'var(--danger-color)' : popData.revel_score >= 0.5 ? 'var(--warning-color)' : 'var(--success-color)';
            html += `
                <div class="prediction-item" style="border-left-color: ${revelColor};">
                    <h4>REVEL</h4>
                    <div class="prediction-score" style="color: ${revelColor};">${popData.revel_score.toFixed(3)}</div>
                    <div><strong>Predizione:</strong> ${popData.revel_score >= 0.75 ? 'Pathogenic' : popData.revel_score >= 0.5 ? 'Likely pathogenic' : 'Uncertain/Benign'}</div>
                    <div style="margin-top: 8px; font-size: 0.9rem;">>0.5 = likely pathogenic</div>
                    <div style="margin-top: 10px;"><a href="https://sites.google.com/site/revelgenomics/" target="_blank" style="color: var(--primary-color);">Vai al tool →</a></div>
                </div>
            `;
        }

        // VEST4
        if (popData.vest4_score !== null && popData.vest4_score !== undefined) {
            html += `
                <div class="prediction-item">
                    <h4>VEST4</h4>
                    <div class="prediction-score">${popData.vest4_score.toFixed(3)}</div>
                    <div style="margin-top: 8px; font-size: 0.9rem;">Functional impact score</div>
                </div>
            `;
        }

        html += '</div>';

        // Add gnomAD frequency section
        if (popData.gnomad_genome_af !== null || popData.gnomad_exome_af !== null) {
            html += '<h3 style="margin-top: 25px;">📊 Frequenze Popolazione (gnomAD)</h3>';
            html += '<div class="prediction-grid">';

            if (popData.gnomad_genome_af !== null && popData.gnomad_genome_af !== undefined) {
                const af = popData.gnomad_genome_af;
                const afFormatted = af < 0.0001 ? af.toExponential(2) : af.toFixed(6);
                html += `
                    <div class="prediction-item">
                        <h4>gnomAD Genomes</h4>
                        <div class="prediction-score">${afFormatted}</div>
                        <div style="margin-top: 8px; font-size: 0.9rem;">Allele Frequency</div>
                    </div>
                `;
            }

            if (popData.gnomad_exome_af !== null && popData.gnomad_exome_af !== undefined) {
                const af = popData.gnomad_exome_af;
                const afFormatted = af < 0.0001 ? af.toExponential(2) : af.toFixed(6);
                html += `
                    <div class="prediction-item">
                        <h4>gnomAD Exomes</h4>
                        <div class="prediction-score">${afFormatted}</div>
                        <div style="margin-top: 8px; font-size: 0.9rem;">Allele Frequency</div>
                    </div>
                `;
            }

            html += `
                <div class="prediction-item" style="grid-column: 1 / -1; background: var(--background); border-left: 4px solid var(--primary-color);">
                    <p style="margin: 0; font-size: 0.95rem;">
                        <strong>Fonte dati:</strong> ${popData.data_source || 'MyVariant.info'}<br>
                        <strong>Note:</strong> Frequenze alleliche da gnomAD v3. AF > 0.05 suggerisce variante comune (BA1). AF > 0.01 suggerisce variante polimorfica (BS1).
                    </p>
                </div>
            `;

            html += '</div>';
        }
    }

    // Missense predictions
    if (predictions.missense && !predictions.missense.error) {
        html += '<h3>Predittori Missense</h3><div class="prediction-grid">';

        for (const [tool, result] of Object.entries(predictions.missense)) {
            if (tool === 'aggregate') continue;
            if (result.error) continue;

            html += `<div class="prediction-item">`;
            html += `<h4>${result.tool || tool}</h4>`;

            if (result.score !== null && result.score !== undefined) {
                html += `<div class="prediction-score">${result.score.toFixed(3)}</div>`;
            }
            if (result.phred_score !== null && result.phred_score !== undefined) {
                html += `<div class="prediction-score">PHRED: ${result.phred_score.toFixed(1)}</div>`;
            }
            if (result.prediction) {
                html += `<div><strong>Predizione:</strong> ${result.prediction}</div>`;
            }
            if (result.interpretation) {
                html += `<div style="margin-top: 8px; font-size: 0.9rem;">${result.interpretation}</div>`;
            }
            if (result.url) {
                html += `<div style="margin-top: 10px;"><a href="${result.url}" target="_blank" style="color: var(--primary-color);">Vai al tool →</a></div>`;
            }

            html += '</div>';
        }

        // Add aggregate
        if (predictions.missense.aggregate) {
            const agg = predictions.missense.aggregate;
            html += `<div class="prediction-item" style="grid-column: 1 / -1; border-left-color: var(--success-color);">`;
            html += `<h4>Consenso</h4>`;
            html += `<div><strong>${agg.consensus}</strong></div>`;
            html += `<div style="margin-top: 8px;">${agg.pathogenic_predictors}/${agg.total_predictors} predittori suggeriscono patogenicità</div>`;
            if (agg.details && agg.details.length > 0) {
                html += '<ul style="margin-top: 10px;">';
                agg.details.forEach(detail => {
                    html += `<li>${detail}</li>`;
                });
                html += '</ul>';
            }
            html += '</div>';
        }

        html += '</div>';
    }

    // Splicing predictions
    if (predictions.splicing && !predictions.splicing.error) {
        html += '<h3 style="margin-top: 25px;">Predittori Splicing</h3><div class="prediction-grid">';

        for (const [tool, result] of Object.entries(predictions.splicing)) {
            if (tool === 'aggregate') continue;
            if (result.error) continue;

            html += `<div class="prediction-item">`;
            html += `<h4>${result.tool || tool}</h4>`;

            if (result.max_delta !== null && result.max_delta !== undefined) {
                html += `<div class="prediction-score">Δ = ${result.max_delta.toFixed(3)}</div>`;
            }
            if (result.interpretation) {
                html += `<div style="margin-top: 8px;">${result.interpretation}</div>`;
            }
            if (result.url) {
                html += `<div style="margin-top: 10px;"><a href="${result.url}" target="_blank" style="color: var(--primary-color);">Vai al tool →</a></div>`;
            }

            html += '</div>';
        }

        // Add aggregate
        if (predictions.splicing.aggregate) {
            const agg = predictions.splicing.aggregate;
            html += `<div class="prediction-item" style="grid-column: 1 / -1; border-left-color: var(--success-color);">`;
            html += `<h4>Consenso Splicing</h4>`;
            html += `<div><strong>${agg.consensus}</strong></div>`;
            if (agg.recommendation) {
                html += `<div style="margin-top: 10px; font-style: italic;">${agg.recommendation}</div>`;
            }
            html += '</div>';
        }

        html += '</div>';
    }

    // Synonymous predictions
    if (predictions.synonymous && !predictions.synonymous.error) {
        html += '<h3 style="margin-top: 25px;">Predittori Varianti Sinonime</h3><div class="prediction-grid">';

        for (const [tool, result] of Object.entries(predictions.synonymous)) {
            if (tool === 'aggregate') continue;
            if (result.error) continue;

            html += `<div class="prediction-item">`;
            html += `<h4>${result.tool || tool}</h4>`;

            if (result.url) {
                html += `<div style="margin-top: 4px;"><a href="${result.url}" target="_blank" style="color: var(--primary-color); font-size: 0.85rem;">${result.url.replace('https://', '')}</a></div>`;
            }
            if (result.interpretation) {
                html += `<div style="margin-top: 8px; font-size: 0.9rem;">${result.interpretation}</div>`;
            }
            if (result.classification) {
                html += `<div style="margin-top: 6px;"><strong>Classificazione MobiDetails:</strong> ${result.classification}</div>`;
            }
            if (result.note) {
                html += `<div style="margin-top: 6px; color: var(--text-secondary); font-size: 0.85rem;">${result.note}</div>`;
            }

            // ESEfinder extra links
            if (result.hexplorer_url) {
                html += `<div style="margin-top: 8px;"><a href="${result.hexplorer_url}" target="_blank" style="color: var(--primary-color);">HEXplorer →</a></div>`;
            }

            html += '</div>';
        }

        // Synonymous aggregate
        if (predictions.synonymous.aggregate) {
            const agg = predictions.synonymous.aggregate;
            const isSyn = agg.is_synonymous;
            const borderColor = isSyn ? 'var(--warning-color)' : 'var(--success-color)';
            html += `<div class="prediction-item" style="grid-column: 1 / -1; border-left-color: ${borderColor};">`;
            html += `<h4>Consenso Variante Sinonima</h4>`;
            html += `<div><strong>${agg.consensus}</strong></div>`;
            if (agg.acmg_criteria && agg.acmg_criteria.length > 0) {
                html += `<div style="margin-top: 8px;"><strong>Criteri ACMG:</strong> ${agg.acmg_criteria.join(', ')}</div>`;
            }
            if (agg.recommendation) {
                html += `<div style="margin-top: 10px; font-style: italic; font-size: 0.9rem;">${agg.recommendation}</div>`;
            }
            html += '</div>';
        }

        html += '</div>';
    }

    // UTR predictions
    if (predictions.utr && !predictions.utr.error) {
        const utrRegion = predictions.utr.aggregate?.utr_region || 'UTR';
        html += `<h3 style="margin-top: 25px;">Predittori Varianti ${utrRegion}</h3><div class="prediction-grid">`;

        // RegulomeDB
        const regulome = predictions.utr.RegulomeDB;
        if (regulome && !regulome.error) {
            const scoreColor = regulome.score && (regulome.score.startsWith('1') || regulome.score.startsWith('2'))
                ? 'var(--danger-color)'
                : regulome.score && (regulome.score === '5' || regulome.score === '6' || regulome.score === '7')
                ? 'var(--success-color)'
                : 'var(--warning-color)';

            html += `<div class="prediction-item" style="border-left-color: ${scoreColor};">`;
            html += `<h4>RegulomeDB</h4>`;
            if (regulome.score) {
                html += `<div class="prediction-score" style="color: ${scoreColor};">Score: ${regulome.score}</div>`;
                html += `<div style="margin-top: 4px; font-size: 0.9rem;">${regulome.score_category || ''}</div>`;
            }
            if (regulome.interpretation) {
                html += `<div style="margin-top: 8px; font-size: 0.9rem;">${regulome.interpretation}</div>`;
            }
            if (regulome.url) {
                html += `<div style="margin-top: 10px;"><a href="${regulome.url}" target="_blank" style="color: var(--primary-color);">Vai a RegulomeDB →</a></div>`;
            }
            html += '</div>';
        }

        // UTR5 or UTR3 specific
        const utr5 = predictions.utr.UTR5Analyzer;
        if (utr5 && !utr5.error) {
            html += `<div class="prediction-item">`;
            html += `<h4>5'UTR Analyzer</h4>`;
            const uorf = utr5.uorf_impact || {};
            if (uorf.potential_uorf_creation) {
                html += `<div style="color: var(--danger-color); font-weight: 600;">⚠ Possibile creazione uORF</div>`;
            } else if (uorf.potential_uorf_disruption) {
                html += `<div style="color: var(--warning-color); font-weight: 600;">⚠ Possibile disruzione uORF</div>`;
            }
            if (utr5.interpretation) {
                html += `<div style="margin-top: 8px; font-size: 0.9rem;">${utr5.interpretation}</div>`;
            }
            if (utr5.url) {
                html += `<div style="margin-top: 10px;"><a href="${utr5.url}" target="_blank" style="color: var(--primary-color);">UTRAnnotator (GitHub) →</a></div>`;
            }
            html += '</div>';
        }

        const utr3 = predictions.utr.UTR3Analyzer;
        if (utr3 && !utr3.error) {
            html += `<div class="prediction-item">`;
            html += `<h4>3'UTR Analyzer</h4>`;
            const polya = utr3.polyadenylation_impact || {};
            if (polya.potential_polya_disruption) {
                html += `<div style="color: var(--danger-color); font-weight: 600;">⚠ Possibile disruzione segnale pA</div>`;
            }
            if (utr3.interpretation) {
                html += `<div style="margin-top: 8px; font-size: 0.9rem;">${utr3.interpretation}</div>`;
            }
            // External links grid
            html += `<div style="margin-top: 10px; display: flex; gap: 8px; flex-wrap: wrap;">`;
            if (utr3.targetscan_url) html += `<a href="${utr3.targetscan_url}" target="_blank" style="color: var(--primary-color); font-size: 0.85rem;">TargetScan →</a>`;
            if (utr3.mirdb_url) html += `<a href="${utr3.mirdb_url}" target="_blank" style="color: var(--primary-color); font-size: 0.85rem;">miRDB →</a>`;
            if (utr3.rbpmap_url) html += `<a href="${utr3.rbpmap_url}" target="_blank" style="color: var(--primary-color); font-size: 0.85rem;">RBPmap →</a>`;
            html += `</div>`;
            html += '</div>';
        }

        // UTR aggregate
        if (predictions.utr.aggregate) {
            const agg = predictions.utr.aggregate;
            html += `<div class="prediction-item" style="grid-column: 1 / -1; border-left-color: var(--primary-color);">`;
            html += `<h4>Consenso ${utrRegion}</h4>`;
            html += `<div><strong>${agg.consensus}</strong></div>`;
            if (agg.acmg_criteria && agg.acmg_criteria.length > 0) {
                html += `<div style="margin-top: 8px;"><strong>Criteri ACMG:</strong> ${agg.acmg_criteria.join(', ')}</div>`;
            }
            if (agg.recommendation) {
                html += `<div style="margin-top: 10px; font-style: italic; font-size: 0.9rem;">${agg.recommendation}</div>`;
            }
            html += '</div>';
        }

        html += '</div>';
    }

    // Indel predictions
    if (predictions.indel && !predictions.indel.error) {
        html += '<h3 style="margin-top: 25px;">Analisi Indel</h3><div class="prediction-grid">';

        const indelAnalysis = predictions.indel.IndelAnalyzer;
        if (indelAnalysis && !indelAnalysis.error) {
            html += `<div class="prediction-item">`;
            html += `<h4>Classificazione Indel</h4>`;
            html += `<div><strong>Tipo:</strong> ${indelAnalysis.variant_type}</div>`;
            html += `<div><strong>Frameshift:</strong> ${indelAnalysis.is_frameshift ? 'Sì' : 'No'}</div>`;
            html += `<div style="margin-top: 10px;">${indelAnalysis.interpretation}</div>`;
            html += '</div>';
        }

        // Add aggregate
        if (predictions.indel.aggregate) {
            const agg = predictions.indel.aggregate;
            html += `<div class="prediction-item" style="border-left-color: var(--success-color);">`;
            html += `<h4>Consenso Indel</h4>`;
            html += `<div><strong>${agg.consensus}</strong></div>`;
            if (agg.acmg_criteria && agg.acmg_criteria.length > 0) {
                html += `<div style="margin-top: 8px;"><strong>Criteri ACMG:</strong> ${agg.acmg_criteria.join(', ')}</div>`;
            }
            if (agg.recommendation) {
                html += `<div style="margin-top: 10px; font-style: italic;">${agg.recommendation}</div>`;
            }
            html += '</div>';
        }

        html += '</div>';
    }

    container.innerHTML = html || '<p>Nessuna predizione disponibile.</p>';
}

function displayLiterature(literature) {
    const container = document.getElementById('literature');
    let html = '';

    if (literature.error) {
        html = `<p>Errore nella ricerca letteratura: ${literature.error}</p>`;
    } else {
        html += `<p><strong>Totale articoli trovati:</strong> ${literature.total_articles}</p>`;

        if (literature.search_url) {
            html += `<p><a href="${literature.search_url}" target="_blank" class="link-button" style="display: inline-block; width: auto;">Cerca su PubMed →</a></p>`;
        }

        // Functional studies
        if (literature.functional_studies && literature.functional_studies.length > 0) {
            html += '<h3 style="margin-top: 25px;">Studi Funzionali</h3>';
            literature.functional_studies.forEach(article => {
                html += formatArticle(article, 'functional');
            });
        }

        // Case reports
        if (literature.case_reports && literature.case_reports.length > 0) {
            html += '<h3 style="margin-top: 25px;">Case Report</h3>';
            literature.case_reports.forEach(article => {
                html += formatArticle(article, 'case');
            });
        }

        // Reviews
        if (literature.reviews && literature.reviews.length > 0) {
            html += '<h3 style="margin-top: 25px;">Review</h3>';
            literature.reviews.forEach(article => {
                html += formatArticle(article, 'review');
            });
        }
    }

    container.innerHTML = html || '<p>Nessun articolo trovato.</p>';
}

function formatArticle(article, type) {
    const badgeClass = type === 'functional' ? 'evidence-functional' :
                       type === 'case' ? 'evidence-case' : 'evidence-review';

    return `
        <div class="literature-item">
            <span class="evidence-badge ${badgeClass}">${type.toUpperCase()}</span>
            <h4>${article.title}</h4>
            <div class="literature-meta">
                ${article.authors} | ${article.journal} | ${article.pub_date}
            </div>
            ${article.abstract ? `<div class="literature-abstract">${article.abstract.substring(0, 300)}...</div>` : ''}
            <a href="${article.pubmed_url}" target="_blank" style="color: var(--primary-color); font-weight: 600;">
                PMID: ${article.pmid} →
            </a>
        </div>
    `;
}

function displayACMG(acmg) {
    const container = document.getElementById('acmgClassification');
    let html = '';

    html += `<div style="background: var(--background); padding: 15px; border-radius: 8px; margin-bottom: 20px;">`;
    html += `<h3>Classificazione: ${acmg.classification}</h3>`;
    html += '</div>';

    // Pathogenic evidence
    if (acmg.pathogenic && acmg.pathogenic.length > 0) {
        html += '<div class="acmg-section">';
        html += '<h4 style="color: var(--danger-color);">Evidenze Patogeniche</h4>';
        acmg.pathogenic.forEach(ev => {
            html += `
                <div class="acmg-item">
                    <span class="acmg-code">${ev.code}</span>
                    <span class="acmg-strength">${ev.strength}</span>
                    <div style="margin-top: 8px;">${ev.description}</div>
                    ${ev.studies ? formatStudies(ev.studies) : ''}
                </div>
            `;
        });
        html += '</div>';
    }

    // Benign evidence
    if (acmg.benign && acmg.benign.length > 0) {
        html += '<div class="acmg-section">';
        html += '<h4 style="color: var(--success-color);">Evidenze Benigne</h4>';
        acmg.benign.forEach(ev => {
            html += `
                <div class="acmg-item">
                    <span class="acmg-code">${ev.code}</span>
                    <span class="acmg-strength">${ev.strength}</span>
                    <div style="margin-top: 8px;">${ev.description}</div>
                </div>
            `;
        });
        html += '</div>';
    }

    if (acmg.pathogenic.length === 0 && acmg.benign.length === 0) {
        html += '<p>Nessuna evidenza ACMG disponibile al momento. I dati di popolazione e gli studi funzionali sono necessari per una classificazione completa.</p>';
    }

    container.innerHTML = html;
}

function formatStudies(studies) {
    if (!studies || studies.length === 0) return '';

    let html = '<div style="margin-top: 10px; padding-left: 20px;">';
    html += '<strong>Studi rilevanti:</strong><ul style="margin-top: 5px;">';
    studies.forEach(study => {
        html += `<li><a href="${study.url}" target="_blank">${study.title} (PMID: ${study.pmid})</a></li>`;
    });
    html += '</ul></div>';
    return html;
}

function displayExternalLinks(links) {
    const container = document.getElementById('externalLinks');

    let html = '<div class="links-grid">';
    for (const [name, url] of Object.entries(links)) {
        const displayName = name.replace(/_/g, ' ');
        html += `<a href="${url}" target="_blank" class="link-button">${displayName}</a>`;
    }
    html += '</div>';

    container.innerHTML = html;
}
