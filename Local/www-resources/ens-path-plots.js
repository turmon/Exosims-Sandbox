/*
  Javascript base for ensemble path plots
  Uses plot.ly JS class to make the plots
  Loads CSV and allows selection of point shading based on its contents.
  
  turmon sep 2018

  TODO:
  (lots)

*/

// this is expected to be defined by a <script> tag in the HTML
if (typeof ens_path_root === "undefined" ) {
   var ens_path_root = '../path-ens';
}
// mode is 'ensemble' or 'single', affects title/legends
if (typeof ens_path_mode === "undefined" ) {
   var ens_path_mode = 'ensemble';
}
//console.log(ens_path_root);

//var sim_url = '//localhost:8080/reduce-ex/reduce-info.csv';
var slews_url = ens_path_root + '/path-slews.csv'
var visits_url = ens_path_root + '/path-visits.csv'

// list of alternate names for the nearest stars
// used only as a display aid
var starAlternateNames2 = {
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

/*
 * A few basic map-based utilities, to enable plotting points midway 
 * along a great circle.  This allows hover-text to be put onto the
 * slews.  Plot.ly does not allow hover-text on lines in maps, so 
 * putting hover-text onto a marker partway along the great circle
 * path is a work-around.
 */
function toDegrees (angle) {
  return angle * (180.0 / Math.PI);
}
function toRadians (angle) {
  return angle * (Math.PI / 180.0);
}
function lonlat2xyz(lon, lat) {
    return {
	x: Math.cos(toRadians(lat)) * Math.cos(toRadians(lon)),
	y: Math.cos(toRadians(lat)) * Math.sin(toRadians(lon)),
	z: Math.sin(toRadians(lat))
    };
}
function xyz2lonlat(pt) {
    var r = Math.sqrt(pt.x*pt.x + pt.y*pt.y + pt.z*pt.z);
    var lon = (toDegrees(Math.atan2(pt.y, pt.x)) + 360) % 360
    var lat = toDegrees(Math.asin(pt.z / r));
    return {lon, lat}
}
function xyz_norm(pt) {
    return Math.sqrt(pt.x*pt.x + pt.y*pt.y + pt.z*pt.z)
}
function gcp1(lon1, lat1, lon2, lat2) {
    // GCP1 = one great circle point in the line from (lon1,lat1) -> (lon2,lat2).
    var p1 = lonlat2xyz(lon1, lat1);
    var p2 = lonlat2xyz(lon2, lat2);
    var delta = {x: p2.x-p1.x, y: p2.y-p1.y, z: p2.z-p1.z};
    // length of the line: between 0 and 2, generally < 1
    var r = xyz_norm(delta);
    // the constant here determines how close the marker appears to the slew endpoint,
    // relative to the total great-circle path length - it is empirical, and the
    // sqrt(r) makes it work for long and short arcs.
    var d = 0.15 / Math.sqrt(r);
    // place one point along the line connecting p1 and p2
    var p_mid = {x:p2.x-d*delta.x, y:p2.y-d*delta.y, z:p2.z-d*delta.z};
    // convert this point to lon,lat
    var {lon, lat} = xyz2lonlat(p_mid);
    // returns an indicator of the length of the great circle, and the
    // lon and lat of a point partway along.
    return {r, lon, lat};
}

// We need two CSV's, but this is not the most elegant way to do it.
Plotly.d3.csv(visits_url, function(err1, visitRows) {
    if (err1) { console.error('Error loading visits CSV: ' + err1); }
    Plotly.d3.csv(slews_url, function(err2, slewRows) { 
	if (err2) { console.error('Error loading slews CSV: ' + err2); }

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
			unit: '[unknown]',
			xform: passthru})
		}
	    }
	    return qoi_info
	}
	    
	// plot the slews and visits
	plotFromSlews('slewPlotDiv', visitRows, slewRows);

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
	function plotFromSlews(target, visitRows, slewRows) {
	    // Make a list of all slews
	    var slew_max = slewRows.reduce(function(a,b) {return Math.max(a,parseFloat(b.slews))}, 0.0);
	    var thick_max = 16.0; // maximum thickness of displayed lines
	    var waypoints = [];
	    var slew_traces = []
	    // primitive equator - no real ticks/tick legends on map objects
	    slew_traces.push({
		lon: [0, 90, 180, 270, 360],
		lat: [0,  0,   0,   0,   0],
		type: 'scattergeo',
		hoverinfo: 'none',
		mode: 'lines',
		line: { color: 'rgba(0, 0, 0, 1)', dash: 'dashdot', width: 0.5 }
	    });
	    // Accumulate the slews (line segments) and waypoints (points)
	    for (var i = 0; i < slewRows.length; i++) {
		var row = slewRows[i];
		var slews = parseFloat(row.slews)
		var label = row.label || ('Slews: ' + slews)
		if (slews < 1.0) {
		    continue; // alter threshold to simplify the plot (disabled)
		}
		//console.log(row.source, row.dest)
		var res = gcp1(visitRows[row.source].lon, visitRows[row.source].lat,
			       visitRows[row.dest].lon,   visitRows[row.dest].lat);
		// we can only show relatively long slews, with labels
		if (res.r > 0.15 && slews/slew_max > 0.08 && label) {
		    res.text = [label,
				'From: ' + visitRows[row.source].name,
				'To: ' + visitRows[row.dest].name].join('<br>');
		    waypoints.push(res)
		}
		slew_traces.push({
		    lon: [visitRows[row.source].lon, visitRows[row.dest].lon],
		    lat: [visitRows[row.source].lat, visitRows[row.dest].lat],
		    type: 'scattergeo',
		    // Note: plotly displays the tooltip at the endpoint, not on the
		    // line itself.  The endpoint text will conflict with any text
		    // on the star itself, so we do not attempt to display it here.
		    hoverinfo: 'none', // suppress the tooltip
		    // text: slews + ' slews from ' + visitRows[row.source].name,
		    mode: 'lines',
		    opacity: 0.4, // there are lines atop lines
		    line: {
			// NB: line color cannot be controlled via a colormap in plotly.
			// You can set it directly - we have not done this because the
			// alternating colors are easier to distinguish.
			//color: 'rgba(55, 128, 191, 0.7)',
			width: (slews/slew_max)*thick_max,
		    },
		})
	    }
	    // Make (x, y, label) for each waypoint
	    var x_way = [], y_way = [], text_way = [];
	    for (var i = 0; i < waypoints.length; i++) {
		x_way.push(waypoints[i].lon);
		y_way.push(waypoints[i].lat);
		text_way.push(waypoints[i].text);
	    }
	    // Make (x, y, visits) for each star
	    var x = [], y = [], z = [], size = [], text = [];
	    for (var i = 0; i < visitRows.length; i++) {
		var row = visitRows[i]
		var visits = parseFloat(row.visit_mean)
		var v_std  = parseFloat(row.visit_std)
		var lon = parseFloat(row.lon)
		var lat = parseFloat(row.lat)
		x.push(lon)
		y.push(lat)
		z.push(visits === 0 ? NaN : visits)
		size.push((visits === 0) ? 3 : 8)
		text.push(row.name + formatStarAlternates(starAlternateNames2, row.name) + '<br>' +
			  visits.toFixed(3) + ' visits' + '<br>' +
			  lon.toFixed(1) + '&deg;, ' + lat.toFixed(1) + '&deg;' )
		//text.push(row['star_Name'] + formatStarAlternates(starAlternateNames, row['star_Name']) + '<br>' +
		//	  qoi_info.name + ' = ' + qoi_value.toPrecision(4) + ' ' + qoi_info.unit + '<br>')
	    }
	    insertPlotly(target, slew_traces,
			 x, y, z, size, text,
			 x_way, y_way, text_way);
	}

	/* Make and insert a plot.ly object into the container in the HTML.
	 * Package the data into JS objects that are passed to plot.ly.
	 * Inputs:
	 *   target -- target HTML container
	 *   slew_traces -- line segments, already as a plot.ly object
	 *   x,y,z,size,text -- star plot attributes, placed into a plot.ly object here
	 *   x_way, y_way, text_way -- waypoint attributes, placed into a plot.ly object here
	 * TODO: Abstract this better.  One way would be to pass in a few objects,
	 * each representing a plot type, rather than separate x,y,z, etc., variables.
	 */
	function insertPlotly(target, slew_traces, x, y, z, size, text, x_way, y_way, text_way) {
	    // var plotDiv = document.getElementById('plotDiv');

	    // default title labels
	    var legend_cbar = 'Mean Cumulative Characterizations';
	    var plot_title = [
		'Mean Cumulative Characterizations and Slew Paths Over Ensemble',
		'Point Shading: Mean Visits across Ensemble; Slew Path Shading: Arbitrary',
		'Diamond Indicators Show Slew Information near Slew Destination'
	    ].join('<br>');
	    if (ens_path_mode === 'single') {
		legend_cbar = 'Number of Characterizations';
		plot_title = [
		    'All Attempted Characterizations and Slew Paths',
		    'Point Shading: Number of Visits; Slew Path Shading: Arbitrary',
		    'Diamond Indicators Show Slew Information near Slew Destination'
		].join('<br>');
	    } 

	    // holds the waypoints
	    var waypoints = {
		type: 'scattergeo',
		lon: x_way,
		lat: y_way,
		text: text_way,
		mode: 'markers',
		hoverinfo: 'text',
		marker: {
		    color: 'rgba(30,30,30,0.4)',
		    size: 7,
		    // not all advertised symbols are available in scattergeo
		    symbol: 'diamond-open',
		    //opacity: 0.7, // opacity is set in "color" above
		},
	    };
	    if (x_way.length > 0) {
		slew_traces.push(waypoints);
	    }
	    // holds the plotted points (stars)
	    var traces = {
		type: 'scattergeo',
		// type: 'scatter',
		lon: x,
		lat: y,
		text: text,
		mode: 'markers',
		hoverinfo: 'text',
		marker: {
		    comment_ : 'leave cmin/cmax unspecified to auto-range',
		    colorscale: 'Portland',
		    color: z,
		    size: size,
		    opacity: 0.7,
		    colorbar: {
			tickfont: {
			    family: 'Arial',
			    size: 16,
			    color: 'black'
			},
			title: '<br\>' + legend_cbar,
			titleside: 'right',
			titlefont: {
			    size: 16,
			    family: 'Arial Bold'
			}
		    }
		},
	    };
	    slew_traces.push(traces)
	    // geographic projection
	    var geo = {
		//resolution: 50,
		showland: false,
		showcountries: false,
		showlakes: false,
		showocean: false,
		projection: {
		    type: 'equirectangular'
		},
		coastlinewidth: 0,
		lataxis: {
		    range: [-80, 80 ],
		    showgrid: true,
		    // axis titles do not work, alas
		    title: 'Longitude [deg]',
		},
		lonaxis:{
		    range: [0, 360],
		    showgrid: true,
		    title: 'Longitude [deg]',
		}
	    };

	    // holds the overall plot layout
	    var layout = {
		title: plot_title, 
		titlefont: {
		    family: "Arial Black",
		    size: 18,
		    color: 'black'
		},
		showlegend: false, // slew legend at right
		geo: geo,
		hovermode: 'closest',
	    }
	    // make the plot from the above objects
	    Plotly.newPlot(target, slew_traces, layout)
	}
    }) // end inner CSV-load
});
