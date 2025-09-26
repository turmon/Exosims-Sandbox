// shim for tabulator.js tables

// endpoints can fail without fatal error, see below
const endpoints = [
    "reduce-yield-all.json", // char data
    "index-files-byname.json", // file counts
    "s_index.json" 
];

function getJSON(response) {
  if (response.ok) return response.json();

  const { status } = response;
    if (false) {
	return {"error": status.toString(10)};
    } else {
	console.log(status);
	const error = new Error(status.toString(10));
	return Promise.reject(error);
    }
}

async function getData(endpoint) {
  const response = await fetch(endpoint);
  return getJSON(response);
}

const data = endpoints.map(getData);

//Promise.allSettled(data).then(results => {
//  console.log('All results:', results);
//});


// needed to get the bottom-row formatting I sought
var summaryNameCalc = function(values, data, calcParams){
    // values - array of column values
    // data - all table data
    // calcParams - params passed from the column definition object
    return `SUMMARY (${values.length} items)`;
};

// last-by-alpha date: for bottom-row computation
var summaryDateCalc = function(values, data, calcParams){
    // values - array of column values
    // data - all table data
    // calcParams - params passed from the column definition object
    return values.sort().pop();
};

var trialInfoTip = function(cell, formatterParams, onRendered){
    // cell - the cell component
    // formatterParams - parameters set for the column
    // onRendered - function to call when the formatter has been rendered
    
    return "Mr" + cell.getValue(); //return the contents of the cell;
}

// custom min/max header filter
var minMaxFilterEditor = function(cell, onRendered, success, cancel, editorParams){

    var end;

    var container = document.createElement("span");

    // create and style inputs
    var start = document.createElement("input");
    start.setAttribute("type", "number");
    start.setAttribute("placeholder", "Min");
    // start.setAttribute("min", 0);
    // start.setAttribute("max", 100);
    start.setAttribute("step", "any");
    start.style.padding = "4px";
    start.style.width = "50%";
    start.style.boxSizing = "border-box";

    start.value = cell.getValue();

    function buildValues(){
        success({
            start:start.value,
            end:end.value,
        });
    }

    function keypress(e){
        if (e.keyCode == 13) {
            buildValues(); // return
        }
        if (e.keyCode == 27) {
            cancel(); // ESC
        }
    }

    end = start.cloneNode();
    end.setAttribute("placeholder", "Max");

    start.addEventListener("change", buildValues);
    start.addEventListener("blur", buildValues);
    start.addEventListener("keydown", keypress);

    end.addEventListener("change", buildValues);
    end.addEventListener("blur", buildValues);
    end.addEventListener("keydown", keypress);

    container.appendChild(start);
    container.appendChild(end);
    return container;
 }

// custom filter function for numeric columns
function minMaxFilterFunction(headerValue, rowValue, rowData, filterParams){
    // headerValue - the value of the header filter element
    // rowValue - the value of the column in this row
    // rowData - the data for the row being filtered
    // filterParams - params object passed to the headerFilterFuncParams property
    // return a boolean: true iff it passes the filter.

    if (rowValue) {
        if (headerValue.start != "") {
            if (headerValue.end != "") {
		// have min and max: filter on both
                return rowValue >= headerValue.start && rowValue <= headerValue.end;
            } else {
		// min, but no max: filter on min only
                return rowValue >= headerValue.start;
            }
        } else {
            if (headerValue.end != "") {
		// max, but no min: filter on max only
                return rowValue <= headerValue.end;
            }
        }
    }
    return true; 
}


Promise.allSettled(data)
    .then(results => {
	let [resReduce, resFiles, resIndex] = results;
	// index fail is OK
	let dataIndex = [{}];
	if (resIndex.status == "fulfilled") {
	    dataIndex = resIndex.value;
	}

	let dataReduce = resReduce.value;
	let dataFiles = resFiles.value;

	// parameter list (empty if s_index.json not present)
	// separate into numerical versus string params
	let dataIndex0 = dataIndex[0];
	const params = Object.keys(dataIndex0).filter(p => {
	    return !p.endsWith("_name");});
	const params_num = params.filter(p => {
	    return typeof dataIndex0[p] === 'number';});
	const params_str = params.filter(p => {
	    return typeof dataIndex0[p] !== 'number';});

	// Cmd-click hides column
	function handleHeaderClick(e, column) {
            if (e.ctrlKey || e.metaKey) {
                // Ctrl+Click or Cmd+Click to hide column
                e.preventDefault();
		column.hide();
                // toggleColumn(column);
                return false; // Prevent sorting
            }
            // Regular click will trigger sorting (handled by Tabulator)
        }

	// displays for numerical and string parameters
	const param_col_num = params_num.map(p => {
	    return {title:`&theta;: ${p}`, field:p,
		    headerClick: handleHeaderClick,
		    headerWordWrap:true,
		    headerFilter:minMaxFilterEditor,
		    headerFilterFunc:minMaxFilterFunction,
		    headerFilterLiveFilter:false,
		   };});
	const param_col_str = params_str.map(p => {
	    return {title:`&theta;: ${p}`, field:p, 
		    headerClick: handleHeaderClick,
		    headerWordWrap:true,
		    headerFilter:"input",
		   };});

	// Map because may not be a 1:1 correspondence
	const newMap = new Map(
	    dataFiles.map(item => {
		// exclude experiment property from map value
		var { experiment, ...newItem } = item;
		// key to announce presence of README information
		if (item.readme_info) {
		    newItem.readmeTag = "&#9432;";
		} else {
		    newItem.readmeTag = "";
		}
		return [experiment, newItem];
	    }
			 ))
		       
	// skip the info column is there's nothing present
	const hasSomeReadmeInfo = dataFiles.some(item => (
	    item.hasOwnProperty("readme_info") && (item.readme_info !== "")))

	// console.log(`Found README info: ${hasSomeReadmeInfo}`)

	// augment dataReduce with dataFiles, as possible
	// NOTE: NaNs in .jsons will cause parse errors upon decoding
	// TODO: handle error if not paired
	
	const dataPlus = dataReduce.map(item => {
	    const pairedItem = newMap.get(item.experiment);
	    return {...item, ...pairedItem}
	});

	// for yields - parameters for the "money" formatter
	const yieldFmtParams = {
	    decimal:".",
	    thousand:false,
	    negativeSign:true,
	    precision:3,
	};

	// vertical height for table (+5 for header/footer, ~25pix/line)
	// (allocate at most this space, but less if the table is short)
	const height_alloc = Math.min(750, (dataPlus.length + 5) * 25);

	var table = new Tabulator("#scenario-table", {
	    height:height_alloc, // set height (CSS or here) to enable Virtual DOM
            data: dataPlus, // assign data to table
            // sortOrderReverse: true, // does not seem to work?
	    layout: "fitColumns", // fit columns to width of table
	    columns:[
		{title:"Name", field:"experiment",
		 width: 350,
		 formatter:"link",
		 formatterParams:{
		     labelField:"experiment",
		     urlField:"url",
		 },
		 headerHozAlign:"center",
		 bottomCalc:summaryNameCalc,
		 headerFilter:"input"},
		// hidden info for readme
		{title:"HiddenInfo", field:"readme_info", visible:false},
		// handle for README
		{title:"&#9432;", field:"readmeTag",
		 visible: hasSomeReadmeInfo,
		 headerSort:false,
		 headerHozAlign:"center",
		 hozAlign:"center",
		 formatter:"html",
		 width:15,
		 clickPopup:function(e, c, onR){return "&#9432;: " + c.getRow().getData().readme_info},
		},
		{title:"Ens. Size", field:"ensemble_size",
		 headerWordWrap:true,
		 headerClick: handleHeaderClick,
		 bottomCalc:"sum"},
		{title:"Last Sim. Date", field:"simtime", hozAlign:"center",
		 headerWordWrap:true,
		 bottomCalc:summaryDateCalc,
		 headerClick: handleHeaderClick,
		 headerFilter:"input"},
		{title:"Reduction Date", field:"runtime", hozAlign:"center",
		 headerWordWrap:true,
		 bottomCalc:summaryDateCalc,
		 headerClick: handleHeaderClick,
		 headerFilter:"input"},
		{title:"User", field:"user",
		 headerClick: handleHeaderClick,
		 headerFilter:"input"},
		{title:"Earths (Det.)", field:"detections_earth_all", hozAlign:"left",
		 headerWordWrap:true,
		 headerClick: handleHeaderClick,
		 headerFilter:minMaxFilterEditor, headerFilterFunc:minMaxFilterFunction,
		 headerFilterLiveFilter:false,
		 formatter:"money", formatterParams:yieldFmtParams,
		 bottomCalcFormatter:"money", bottomCalcFormatterParams:yieldFmtParams,
		 bottomCalc:"max"},
		{title:"Earths (Char.)", field:"chars_earth_unique", hozAlign:"left",
		 headerWordWrap:true,
		 headerClick: handleHeaderClick,
		 headerFilter:minMaxFilterEditor, headerFilterFunc:minMaxFilterFunction,
		 headerFilterLiveFilter:false,
		 formatter:"money", formatterParams:yieldFmtParams,
		 bottomCalcFormatter:"money", bottomCalcFormatterParams:yieldFmtParams,
		 bottomCalc:"max"},
		{title:"Earths (Strict)", field:"chars_earth_strict", hozAlign:"left",
		 headerWordWrap:true,
		 headerClick: handleHeaderClick,
		 headerFilter:minMaxFilterEditor, headerFilterFunc:minMaxFilterFunction,
		 headerFilterLiveFilter:false,
		 formatter:"money", formatterParams:yieldFmtParams,
		 bottomCalcFormatter:"money", bottomCalcFormatterParams:yieldFmtParams,
		 bottomCalc:"max"},
		{title:"Ens. Graphs", field:"index_gfx_count", bottomCalc:"sum",
		 headerClick: handleHeaderClick,
		 headerWordWrap:true,},
		{title:"Path Summ's", field:"index_path_count", bottomCalc:"sum",
		 headerClick: handleHeaderClick,
		 headerWordWrap:true,},
		{title:"Path Graphs", field:"index_path_gfx", bottomCalc:"sum",
		 headerClick: handleHeaderClick,
		 headerWordWrap:true,},
	    ].concat(param_col_str).concat(param_col_num),
	});
	// below: further failed attempt to change sort direction
	//table.setSort([
	//    {column:"simtime", dir:"desc"}, ]);
	//console.log(dataPlus);
    })
    .catch(error => {
	// this will double-log the fetch errors
	console.error('There was a problem fetching the data:', error);
	throw error;
    });


