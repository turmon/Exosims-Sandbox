/*
  Javascript base for radius/luminosity "shaded by X" plots
  Uses plot.ly JS class to make the plots
  Loads two CSV files and allows selection of point shading based on fields within.
  
  This version inserts detection-related plots into one HTML <div>,
  and characterization plots into another one.

  turmon sep 2018

*/


var sim_url = '../reduce-info.csv';
var data_url = '../reduce-star-target.csv';
//var sim_url = '//localhost:8080/reduce/reduce-info.csv';
//var data_url = '//localhost:8080/reduce-ex/reduce-star-target.csv';

function passthru(x) {return x};

// metadata about a selected list of QOIs
//   if we add a new one to the CSV, it might not appear here right away,
//   so not all plot-able QOIs are required to have an entry below
//   screen: fieldname that will treat a given row in the present field as if it's NaN, if the corresponding row in "screen" is NaN
var known_qoi_info = [
    // detection non-yield
    {fieldname: 'h_star_det_visit_mean',        name: 'Detection Visits',                  screen: 'h_star_det_plan_value_mean', unit: 'count', xform: passthru},
    {fieldname: 'h_star_det_comp_mean',         name: 'Mean Completeness (Det.)',          screen: '',                           unit: 'count', xform: passthru},
    {fieldname: 'h_star_det_tInt_mean',         name: 'Mean Cume Integ. Time (Det.)',      screen: 'h_star_det_plan_value_mean', unit: 'day',   xform: Math.log10},
    {fieldname: 'h_star_det_tIntAvg_mean',      name: 'Mean One-visit Integ. Time (Det.)', screen: 'h_star_det_plan_value_mean', unit: 'day',   xform: Math.log10},
    {fieldname: 'h_star_det_tobs1_mean',        name: 'First Observation Time (Det.)',     screen: 'h_star_det_plan_value_mean', unit: 'day',   xform: passthru},
    // detection yields			        
    {fieldname: 'h_star_det_plan_cume_mean',    name: 'Mean Total Detections',             screen: 'h_star_det_plan_value_mean', unit: 'count', xform: passthru},
    {fieldname: 'h_star_det_plan_uniq_mean',    name: 'Mean Unique Detections',            screen: 'h_star_det_plan_value_mean', unit: 'count', xform: passthru},
    {fieldname: 'h_star_det_earth_cume_mean',   name: 'Mean Total Detections: Earths',     screen: 'h_star_det_plan_value_mean', unit: 'count', xform: passthru},
    {fieldname: 'h_star_det_earth_uniq_mean',   name: 'Mean Unique Detections: Earths',    screen: 'h_star_det_plan_value_mean', unit: 'count', xform: passthru},
    // detection ratios			        
    {fieldname: 'h_star_det_plan_value_mean',   name: 'Detection Rank',                    screen: '',                           unit: 'count/day',   xform: Math.log10},
    {fieldname: 'h_star_det_plan_frac_mean',    name: 'Planets Detected/Planets Present',  screen: 'h_star_det_plan_value_mean', unit: 'count/count', xform: passthru},
    {fieldname: 'h_star_det_earth_value_mean',  name: 'Earth Detection Rank',              screen: '',                           unit: 'count/day',   xform: Math.log10},
    {fieldname: 'h_star_det_earth_frac_mean',   name: 'Earths Detected/Earths Present',    screen: 'h_star_det_plan_value_mean', unit: 'count/count', xform: passthru},
    // char non-yield			        
    {fieldname: 'h_star_char_visit_mean',       name: 'Characterization Visits',           screen: 'h_star_char_plan_value_mean', unit: 'count', xform: passthru},
    {fieldname: 'h_star_char_comp_mean',        name: 'Mean Completeness (Char.)',         screen: '',                            unit: 'count', xform: passthru},
    {fieldname: 'h_star_char_tInt_mean',        name: 'Mean Cume Integ. Time (Char.)',     screen: 'h_star_char_plan_value_mean', unit: 'day',   xform: Math.log10},
    {fieldname: 'h_star_char_tIntAvg_mean',     name: 'Mean One-visit Integ. Time (Char.)',screen: 'h_star_char_plan_value_mean', unit: 'day',   xform: Math.log10},
    {fieldname: 'h_star_char_tobs1_mean',       name: 'First Observation Time (Char.)',    screen: 'h_star_char_plan_value_mean', unit: 'day',   xform: passthru},
    // char yields
    {fieldname: 'h_star_char_plan_cume_mean',   name: 'Mean Total Characterizations',          screen: 'h_star_char_plan_value_mean', unit: 'count', xform: passthru},
    {fieldname: 'h_star_char_plan_uniq_mean',   name: 'Mean Unique Characterizations',         screen: 'h_star_char_plan_value_mean', unit: 'count', xform: passthru},
    {fieldname: 'h_star_char_earth_cume_mean',  name: 'Mean Total Characterizations: Earths',  screen: 'h_star_char_plan_value_mean', unit: 'count', xform: passthru},
    {fieldname: 'h_star_char_earth_uniq_mean',  name: 'Mean Unique Characterizations: Earths', screen: 'h_star_char_plan_value_mean', unit: 'count', xform: passthru},
    // char ratios
    {fieldname: 'h_star_char_plan_value_mean',  name: 'Characterization Rank',                 screen: '',                            unit: 'count/day',   xform: Math.log10},
    {fieldname: 'h_star_char_plan_frac_mean',   name: 'Planets Characterized/Planets Present', screen: 'h_star_char_plan_value_mean', unit: 'count/count', xform: passthru},
    {fieldname: 'h_star_char_earth_value_mean', name: 'Earth Characterization Rank',           screen: '',                            unit: 'count/day',   xform: Math.log10},
    {fieldname: 'h_star_char_earth_frac_mean',  name: 'Earths Characterized/Earths Present',   screen: 'h_star_char_plan_value_mean', unit: 'count/count', xform: passthru},
    // others
    {fieldname: 'h_star_plan_per_star_mean',    name: 'Planets per Star',                      screen: '',                            unit: 'count', xform: passthru},
    {fieldname: 'h_star_earth_per_star_mean',   name: 'Earths per Star',                       screen: '',                            unit: 'count', xform: passthru},
    // promotions
    {fieldname: 'h_star_promo_allplan_mean',    name: 'Promotion Rate',                        screen: 'h_star_char_plan_value_mean', unit: 'count', xform: passthru},
    {fieldname: 'h_star_promo_hzone_mean',      name: 'Promotion Rate: Stars with HZ Planets', screen: 'h_star_char_plan_value_mean', unit: 'count', xform: passthru},
    {fieldname: 'h_star_promo_earth_mean',      name: 'Promotion Rate: Stars with Earths',     screen: 'h_star_char_plan_value_mean', unit: 'count', xform: passthru},
];

// list of alternate names for the nearest stars
// used only as a display aid
var starAlternateNames = {
    'HIP 5336':   ['mu Cas'],
    // automatically found, see "star-alt-names-merged.js"
    "HIP 677":  ["HD 358", "LTT 10039", "alpha And"],
    "HIP 1475":  ["HD 1326", "Gl 15A", "GX Andromedae"],
    "HIP 2021":  ["HD 2151", "Gl 19", "beta Hyi"],
    "HIP 3821":  ["HD 4614", "Gl 34A", "eta Cassiopei A"],
    "HIP 3829":  ["GJ 35", "van Maanen's Star"],
    "HIP 5643":  ["GJ 54.1", "YZ Ceti"],
    "HIP 8102":  ["HD 10700", "Gl 71", "tau Ceti"],
    "HIP 9884":  ["HD 12929", "Gl 84.3", "alf Arietis "],
    "HIP 15510":  ["HD 20794", "Gl 139", "82 Eridani"],
    "HIP 16537":  ["HD 22049", "Gl 144", "epsilon Eridani"],
    "HIP 19849":  ["HD 26965", "Gl 166A", "omicron 2 Eridani"],
    "HIP 21088":  ["Gl 169.1A", "Stein 2051"],
    "HIP 21421":  ["HD 29139", "Gl 171.1A", "Aldebaran"],
    "HIP 22449":  ["HD 30652", "Gl 178", "1 Ori"],
    "HIP 24186":  ["HD 33793", "Gl 191", "Kapteyn's Star"],
    "HIP 24608":  ["HD 34029", "Gl 194A", "Capella AB"],
    "HIP 25878":  ["HD 36395", "Gl 205", "Wolf 1453"],
    "HIP 26857":  ["Gl 213", "Ross 47"],
    "HIP 27913":  ["HD 39587", "Gl 222AB", "chi Ori"],
    "HIP 28360":  ["HD 40183", "beta Aur AB"],
    "HIP 30920":  ["Gl 234A", "Ross 614 A"],
    "HIP 32349":  ["HD 48915", "Gl 244A", "Sirius A"],
    "HIP 34603":  ["GJ 268 A", "QY Aurigae A"],
    "HIP 36208":  ["Gl 273", "Luyten's Star"],
    "HIP 36850":  ["HD 60179", "Gl 278A", "Castor"],
    "HIP 37279":  ["HD 61421", "Gl 280A", "Procyon A"],
    "HIP 37766":  ["Gl 285", "Ross 882"],
    "HIP 37826":  ["HD 62509", "Gl 286", "Pollux"],
    "HIP 42913":  ["HD 74956", "Gl 321.3A", "delta Velorum A"],
    "HIP 49669":  ["HD 87901", "LTT 12716", "Regulus"],
    "HIP 53020":  ["GJ 402", "EE Leonis"],
    "HIP 53767":  ["Gl 408", "Ross 104"],
    "HIP 54035":  ["HD 95735", "Gl 411", "Lalande 21185"],
    "HIP 57367":  ["Gl 440", "WD 1142-645"],
    "HIP 57548":  ["Gl 447", "Ross 128"],
    "HIP 61084":  ["HD 108903", "Gl 470", "Gacrux"],
    "HIP 61317":  ["HD 109358", "Gl 475", "beta CVn"],
    "HIP 62956":  ["HD 112185", "eps UMa"],
    "HIP 67155":  ["HD 119850", "Gl 526", "Wolf 498"],
    "HIP 69673":  ["HD 124897", "Gl 541", "Arcturus"],
    "HIP 70890":  ["Gl 551", "prox Cen"],
    "HIP 71253":  ["Gl 555", "HN Librae"],
    "HIP 71681":  ["HD 128621", "Gl 559B", "alpha Cen B"],
    "HIP 71683":  ["HD 128620", "Gl 559A", "alpha Cen A"],
    "HIP 72659":  ["HD 131156", " GJ 566A", "ksi Bootis A"],
    "HIP 74995":  ["Gl 581", "Wolf 562"],
    "HIP 80824":  ["Gl 628", "Wolf 1061"],
    "HIP 82809":  ["Gl 643", "Wolf 629"],
    "HIP 82817":  ["HD 152751", "Gl 644A", "Wolf 630 A"],
    "HIP 84405":  ["HD 155885", "Gl 663A", "36 Ophiuchi A"],
    "HIP 84478":  ["HD 156026", "Gl 664", "36 Ophiuchi C"],
    "HIP 87937":  ["Gl 699", "Barnard's Star"],
    "HIP 88601":  ["HD 165341", "Gl 702A", "70 Ophiuchi A"],
    "HIP 89937":  ["HD 170153", "Gl 713AB", "chi Dra"],
    "HIP 91262":  ["HD 172167", "Gl 721", "Vega"],
    "HIP 92403":  ["GJ 729", "Ross 154"],
    "HIP 94761":  ["HD 180617", "Gl 752A", "Wolf 1055"],
    "HIP 96100":  ["HD 185144", "Gl 764", "sigma Draconis"],
    "HIP 97649":  ["HD 187642", "Gl 768", "Altair"],
    "HIP 99240":  ["HD 190248", "Gl 780", "delta Pavonis"],
    "HIP 104214": ["HD 201091", "Gl 820A", "61 Cygni A"],
    "HIP 104217": ["HD 201092", "Gl 820B", "61 Cygni B"],
    "HIP 105090": ["HD 202560", "Gl 825", "AX Microscopii"],
    "HIP 106106": ["Gl 829", "Ross 775 A"],
    "HIP 108870": ["HD 209100", "Gl 845", "epsilon Indi A"],
    "HIP 110893": ["HD 239960", "GJ 860 A", "Kruger 60 A"],
    "HIP 112460": ["GJ 873", "EV Lacertae"],
    "HIP 113020": ["Gl 876", "Ross 780"],
    "HIP 113296": ["HD 216899", "Gl 880", "Ross 671"],
    "HIP 113368": ["HD 216956", "Gl 881A", "Fomalhaut"],
    "HIP 114046": ["HD 217987", "Gl 887", "Lacaille 9352"],
    "HIP 116132": ["GJ 896 A", "EQ Pegasi"],
    'HIP 439':    ['Gl 1', 'CD-37 15492'],
    'HIP 1242':   ['GJ 1005 A', 'L 722-22 A'],
    'HIP 29295':  ['Gl 229 A', 'L 668-21 A'],
    'HIP 33226':  ['Gl 251', 'Wolf 294'],
    'HIP 49908':  ['Gl 380', 'Groombridge 1618'],
    'HIP 54211':  ['Gl 412 A', 'Lalande 21258 A'],
    'HIP 57544':  ['Gl 445', 'G254-29'],
    'HIP 73182':  ['Gl 570 B', 'Lalande 27173 B'],
    'HIP 73184':  ['Gl 570 A', 'Lalande 27173 A'],
    'HIP 76074':  ['Gl 588', 'CD-40 9712'],
    'HIP 85523':  ['Gl 674', 'CD-46 11540'],
    'HIP 86162':  ['Gl 687', 'BD+68 946'],
    'HIP 86214':  ['Gl 682', 'CD-44 11909'],
    'HIP 86990':  ['Gl 693', 'L 205-128'],
    'HIP 91768':  ['Gl 725 A', 'Struve 2398 A'],
    'HIP 91772':  ['Gl 725 B', 'Struve 2398 B'],
    'HIP 99461':  ['Gl 783 A', 'HR 7703A'],
    'HIP 106440': ['Gl 832', 'L 354-89'],
    'HIP 117473': ['Gl 908', 'Lalande 46650'],
};


Plotly.d3.csv(sim_url, function(err1, simRow) {
    if (err1) { console.error('Error loading sim CSV: ' + err1); }
    Plotly.d3.csv(data_url, function(err2, allRows) { 
	if (err2) { console.error('Error loading CSV: ' + err2); }

	// match fieldname list to known metadata (units, plot name), where possible,
	// using metadata table at top of file
	function matchFieldnamesToInfo(qoi_fieldnames) {
	    var qoi_info = [];
	    for (var i = 0; i < qoi_fieldnames.length; i++) {
		var qoi_fieldname = qoi_fieldnames[i];
		var qoi_matches = known_qoi_info.filter(function(info) {return info.fieldname === qoi_fieldname})
		if (qoi_matches.length === 1) {
		    qoi_info.push(qoi_matches[0])
		} else {
		    // push a generic info block - units are unknown
		    qoi_info.push({
			fieldname: qoi_fieldname,
			name: qoi_fieldname.replace(/_mean$/, '').replace(/^h_star_/, ''),
			screen: '',
			unit: '[unknown]',
			xform: passthru})
		}
	    }
	    return qoi_info
	}
	    
	// set up the selector determining which QOI we will plot
	function assignPlotVariableOptions(qoi_info, selector) {
	    for (var i = 0; i < qoi_info.length;  i++) {
		var currentOption = document.createElement('option');
		currentOption.text = qoi_info[i].name + ' [' + qoi_info[i].fieldname + ']';
		selector.appendChild(currentOption);
	    }
	}
	
	// overall title of the experiment, obtained from outer CSV
	var simBigTitle = 'Experiment ' + simRow[0].experiment + ', Ensemble Size ' + simRow[0].ensemble_size;

	// the fieldnames are unique and define what CSV columns we can plot
	var qoi_gen_fieldnames = Object.getOwnPropertyNames(allRows[0])
	    .filter(function (name) {return name.match(/h_star_.*_per_star_mean/)});
	var qoi_det_fieldnames = Object.getOwnPropertyNames(allRows[0])
	    .filter(function (name) {return name.match(/h_star_(?:det|promo).*_mean/)})
	    .concat(qoi_gen_fieldnames);
	var qoi_char_fieldnames = Object.getOwnPropertyNames(allRows[0])
	    .filter(function (name) {return name.match(/h_star_char.*_mean/)})
	    .concat(qoi_gen_fieldnames);

	// the best metadata we have
	var qoi_det_info  = matchFieldnamesToInfo(qoi_det_fieldnames);
	var qoi_char_info = matchFieldnamesToInfo(qoi_char_fieldnames);

	// initial plot of some QOI
	plotFromRows('detPlotDiv',  qoi_det_info[0],  allRows);
	plotFromRows('charPlotDiv', qoi_char_info[0], allRows);

	// set up the plotted-variable button, and its callback
	var qoiDetSelector = document.querySelector('.det_qoi_select');
	assignPlotVariableOptions(qoi_det_info, qoiDetSelector);
	qoiDetSelector.addEventListener(
	    'change',
	    function () {
		plotFromRows('detPlotDiv', qoi_det_info[qoiDetSelector.selectedIndex], allRows)
	    },
	    false);
	
	var qoiCharSelector = document.querySelector('.char_qoi_select');
	assignPlotVariableOptions(qoi_char_info, qoiCharSelector);
	qoiCharSelector.addEventListener(
	    'change',
	    function () {
		plotFromRows('charPlotDiv', qoi_char_info[qoiCharSelector.selectedIndex], allRows)
	    },
	    false);
	

	function formatStarAlternates(altNames, star) {
	    // strip trailing A/B/C from the HIP NNNN input name: match on just the HIP ID#
	    var star_alone = star.replace(/ [A-C]$/, '')
	    if (altNames.hasOwnProperty(star_alone)) {
		return ' (' + altNames[star_alone].join(', ') + ')'
	    } else {
		return ''
	    }
	}
	/* filter the CSV lines into the information needed for the plot specified in qoi_info,
	 * then make a plot with this information
	 */
	function plotFromRows(target, qoi_info, lines) {
	    // field names
	    var qoi_fieldname_mean = qoi_info.fieldname
	    var qoi_fieldname_std  = qoi_fieldname_mean.replace(/_mean$/, '_std')
	    var qoi_fieldname_sem  = qoi_fieldname_mean.replace(/_mean$/, '_sem')
	    var qoi_fieldname_nens = qoi_fieldname_mean.replace(/_mean$/, '_nEns')
	    var qoi_screen_name = qoi_info.screen
	    // cons up lists of points, colors, sizes, text annotations
	    var x = [], y = [], qoi = [], size = [], text = [];
	    for (var i = 0; i < lines.length; i++) {
		var row = lines[i];
		var qoi_value = parseFloat(row[qoi_fieldname_mean])
		var qoi_nens  = parseInt(  row[qoi_fieldname_nens])
		var qoi_sem   = parseFloat(row[qoi_fieldname_sem])
		var qoi_std   = parseFloat(row[qoi_fieldname_std])
		var qoi_screen = qoi_screen_name ? parseFloat(row[qoi_screen_name]) : 0.0 // anything non-nan is fine
		x.push(row['star_dist'])
		y.push(row['star_L'])
		qoi.push(qoi_info.xform(qoi_value))
		size.push((isNaN(qoi_value) || isNaN(qoi_screen)) ? 5 : 8)
		text.push(`${row['star_Name']} ${formatStarAlternates(starAlternateNames, row['star_Name'])} [${row['star_Spec']}]<br>` +
			  `${qoi_info.name} = ${qoi_value.toPrecision(4)} ${qoi_info.unit} &plusmn; ${qoi_sem.toPrecision(3)} [N = ${qoi_nens}]<br>` +
			  `Std. Deviation = ${qoi_std.toPrecision(4)} ${qoi_info.unit}`)
	    }
	    // console.log(qoi_info)
	    insertPlotly(target, qoi_info, x, y, qoi, size, text);
	}

	/* Make and insert a plot.ly object into the container in the HTML.
	 * Package the data into JS objects that are passed to plot.ly.
	 */
	function insertPlotly(target, qoi_info, x, y, qoi, size, text) {
	    // var plotDiv = document.getElementById('plotDiv');
	    // if xform(1.0) = 1.0, the QOI was not transformed, else it was log10
	    var xform_label = (qoi_info.xform(1.0) === 1.0) ? '' : 'log<sub>10</sub> '
	    // holds the plotted points
	    var traces = [{
		x: x,
		y: y,
		text: text,
		mode: 'markers',
		type: 'scatter',
		hoverinfo: 'text',
		marker: {
		    comment_ : 'leave cmin/cmax unspecified to auto-range',
		    colorscale: 'Portland',
		    color: qoi,
		    size: size,
		    opacity: 0.7,
		    colorbar: {
			tickfont: {
			    family: 'Arial',
			    size: 16,
			    color: 'black'
			},
			title: '<br\>' + qoi_info.name + ' [' + xform_label + qoi_info.unit + ']',
			titleside: 'right',
			titlefont: {
			    size: 16,
			    family: 'Arial Bold'
			}
		    }
		},
	    }];
	    // holds the overall plot layout
	    var layout = {
		title: ['Star Luminosity vs. Distance, ' +
			'Shaded by: ' + qoi_info.name,
			simBigTitle].join('<br>'),
		titlefont: {
		    family: "Arial Black",
		    size: 18,
		    color: 'black'
		},
		xaxis: {
		    title: 'Distance [pc]',
		    titlefont: {
			family: "Arial Bold",
			size: 18,
		    },
		    size: 20,
		    tickfont: {
			family: 'Arial',
			size: 16,
			color: 'black'
		    },
		    linecolor: 'black',
		    linewidth: 2,
		    mirror: true,
		    range: [0.0, 30.0],
		    // rangemode: 'tozero', // go to zero [pc]
		},
		yaxis: {
		    title: 'Luminosity [L<sub>sun</sub>]',
		    titlefont: {
			family: "Arial Bold",
			size: 18,
		    },
		    type: 'log',
		    size: 20,
		    tickfont: {
			family: 'Arial',
			size: 14,
			color: 'black'
		    },
		    linecolor: 'black',
		    linewidth: 2,
		    mirror: true
		},
		hovermode: 'closest',
		autosize: true,
	    }
	    // make the plot from the above objects
	    Plotly.newPlot(target, traces, layout)
	}
    }) // end inner CSV-load
});
