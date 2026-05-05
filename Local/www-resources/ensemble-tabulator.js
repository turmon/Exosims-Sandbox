/*
  Interactive ensemble summary table for an EXOSIMS family-level index page.

  Fetches three JSON files and builds a sortable, filterable Tabulator table
  in #scenario-table:
    reduce-yield-all.json    -- per-ensemble yield statistics (required)
    index-files-byname.json  -- per-ensemble file counts and URLs (required)
    s_index.json             -- per-ensemble parameter values (optional)

  If either required file is missing, a red error message is shown inside the
  table div instead of leaving it silently blank.

  The script tag can be placed in <head> because the entry point is wrapped
  in a DOMContentLoaded listener.

  turmon 2019, revised may 2026
*/

(function () {

const endpoints = [
    'reduce-yield-all.json',    // yield data (required)
    'index-files-byname.json',  // file counts (required)
    's_index.json',             // parameter index (optional)
];

// Write a visible error message into the table div.
function showError(message) {
    const div = document.getElementById('scenario-table');
    if (div)
        div.innerHTML =
            '<p style="color:red; padding:20px; font-size:1.1em;">' +
            '<strong>Could not load table data:</strong> ' + message + '</p>';
}

// Fetch a JSON file; rejects with a descriptive Error on HTTP failure.
function loadJSON(url) {
    return fetch(url).then(function(response) {
        if (!response.ok)
            throw new Error('Could not load ' + url + ' (HTTP ' + response.status + ')');
        return response.json();
    });
}

// Bottom-row summary formatters for the table footer.
function summaryNameCalc(values) {
    return `SUMMARY (${values.length} items)`;
}
function summaryDateCalc(values) {
    return values.sort().pop();
}

// Custom min/max range header filter: builds a {start, end} value object.
function minMaxFilterEditor(cell, onRendered, success, cancel) {
    const container = document.createElement('span');

    const start = document.createElement('input');
    start.setAttribute('type', 'number');
    start.setAttribute('placeholder', 'Min');
    start.setAttribute('step', 'any');
    start.style.padding   = '4px';
    start.style.width     = '50%';
    start.style.boxSizing = 'border-box';
    start.value = cell.getValue();

    const end = start.cloneNode();
    end.setAttribute('placeholder', 'Max');

    function buildValues() {
        success({ start: start.value, end: end.value });
    }
    function keypress(e) {
        if (e.keyCode === 13) buildValues(); // Return
        if (e.keyCode === 27) cancel();      // ESC
    }

    [start, end].forEach(function(inp) {
        inp.addEventListener('change',  buildValues);
        inp.addEventListener('blur',    buildValues);
        inp.addEventListener('keydown', keypress);
    });

    container.appendChild(start);
    container.appendChild(end);
    return container;
}

// Companion filter function: true iff rowValue falls in [start, end].
function minMaxFilterFunction(headerValue, rowValue) {
    if (rowValue != null) {
        if (headerValue.start !== '') {
            if (headerValue.end !== '')
                return rowValue >= headerValue.start && rowValue <= headerValue.end;
            return rowValue >= headerValue.start;
        } else if (headerValue.end !== '') {
            return rowValue <= headerValue.end;
        }
    }
    return true;
}


// Ctrl/Cmd-click on a column header hides that column.
function handleHeaderClick(e, column) {
    if (e.ctrlKey || e.metaKey) {
        e.preventDefault();
        column.hide();
        return false;
    }
}


document.addEventListener('DOMContentLoaded', function() {
    if (!document.getElementById('scenario-table')) return;

    const data = endpoints.map(loadJSON);

    Promise.allSettled(data)
        .then(function(results) {
            const [resReduce, resFiles, resIndex] = results;

            // Both required files must have loaded successfully.
            if (resReduce.status !== 'fulfilled' || resFiles.status !== 'fulfilled') {
                const who = resReduce.status !== 'fulfilled'
                    ? 'reduce-yield-all.json' : 'index-files-byname.json';
                showError('Required file missing: ' + who);
                return;
            }

            const dataReduce = resReduce.value;
            const dataFiles  = resFiles.value;

            // s_index.json is optional; an empty object stands in when absent.
            const dataIndex  = resIndex.status === 'fulfilled' ? resIndex.value : [{}];

            // Split parameter columns into numeric vs. string for appropriate filtering.
            const dataIndex0 = dataIndex[0];
            const params     = Object.keys(dataIndex0).filter(p => !p.endsWith('_name'));
            const params_num = params.filter(p => typeof dataIndex0[p] === 'number');
            const params_str = params.filter(p => typeof dataIndex0[p] !== 'number');

            const param_col_num = params_num.map(p => ({
                title: `&theta;: ${p}`, field: p,
                headerClick: handleHeaderClick,
                headerWordWrap: true,
                headerFilter: minMaxFilterEditor,
                headerFilterFunc: minMaxFilterFunction,
                headerFilterLiveFilter: false,
            }));
            const param_col_str = params_str.map(p => ({
                title: `&theta;: ${p}`, field: p,
                headerClick: handleHeaderClick,
                headerWordWrap: true,
                headerFilter: 'input',
            }));

            // Join file-count data onto yield data by experiment name.
            const fileMap = new Map(
                dataFiles.map(function(item) {
                    const { experiment, ...rest } = item;
                    rest.readmeTag = item.readme_info ? '&#9432;' : '';
                    return [experiment, rest];
                })
            );
            const hasSomeReadmeInfo = dataFiles.some(
                item => item.hasOwnProperty('readme_info') && item.readme_info !== '');

            const dataPlus = dataReduce.map(function(item) {
                const paired = fileMap.get(item.experiment);
                if (!paired) console.warn('ensemble-tabulator: no file data for experiment:', item.experiment);
                return { ...item, ...paired };
            });

            const yieldFmtParams = {
                decimal: '.', thousand: false, negativeSign: true, precision: 3,
            };

            function yieldCol(title, field) {
                return {title, field, hozAlign: 'left', headerWordWrap: true,
                        headerClick: handleHeaderClick,
                        headerFilter: minMaxFilterEditor,
                        headerFilterFunc: minMaxFilterFunction,
                        headerFilterLiveFilter: false,
                        formatter: 'money', formatterParams: yieldFmtParams,
                        bottomCalcFormatter: 'money', bottomCalcFormatterParams: yieldFmtParams,
                        bottomCalc: 'max'};
            }

            // Allocate at most 750px; less if the table is short (~25px per row).
            const height_alloc = Math.min(750, (dataPlus.length + 5) * 25);

            new Tabulator('#scenario-table', {
                height: height_alloc,
                data: dataPlus,
                layout: 'fitColumns',
                columns: [
                    {title: 'Name', field: 'experiment',
                     width: 350,
                     formatter: 'link',
                     formatterParams: { labelField: 'experiment', urlField: 'url' },
                     headerHozAlign: 'center',
                     bottomCalc: summaryNameCalc,
                     headerFilter: 'input'},
                    {title: 'HiddenInfo', field: 'readme_info', visible: false},
                    {title: '&#9432;', field: 'readmeTag',
                     visible: hasSomeReadmeInfo,
                     headerSort: false,
                     headerHozAlign: 'center',
                     hozAlign: 'center',
                     formatter: 'html',
                     width: 15,
                     clickPopup: function(e, c) {
                         return '&#9432;: ' + c.getRow().getData().readme_info;
                     }},
                    {title: 'Ens. Size', field: 'ensemble_size',
                     headerWordWrap: true,
                     headerClick: handleHeaderClick,
                     bottomCalc: 'sum'},
                    {title: 'Last Sim. Date', field: 'simtime', hozAlign: 'center',
                     headerWordWrap: true,
                     bottomCalc: summaryDateCalc,
                     headerClick: handleHeaderClick,
                     headerFilter: 'input'},
                    {title: 'Reduction Date', field: 'runtime', hozAlign: 'center',
                     headerWordWrap: true,
                     bottomCalc: summaryDateCalc,
                     headerClick: handleHeaderClick,
                     headerFilter: 'input'},
                    {title: 'User', field: 'user',
                     headerClick: handleHeaderClick,
                     headerFilter: 'input'},
                    yieldCol('Earths (All)',    'detections_earth_all'),
                    yieldCol('Earths (Det.)',   'detections_earth_unique'),
                    yieldCol('Earths (Char.)',  'chars_earth_unique'),
                    yieldCol('Earths (Strict)', 'chars_earth_strict'),
                    {title: 'Ens. Graphs', field: 'index_gfx_count',
                     headerClick: handleHeaderClick,
                     headerWordWrap: true,
                     bottomCalc: 'sum'},
                    {title: "Path Summ's", field: 'index_path_count',
                     headerClick: handleHeaderClick,
                     headerWordWrap: true,
                     bottomCalc: 'sum'},
                    {title: 'Path Graphs', field: 'index_path_gfx',
                     headerClick: handleHeaderClick,
                     headerWordWrap: true,
                     bottomCalc: 'sum'},
                ].concat(param_col_str).concat(param_col_num),
            });
        })
        .catch(function(error) {
            showError(error.message);
        });
});

}());
