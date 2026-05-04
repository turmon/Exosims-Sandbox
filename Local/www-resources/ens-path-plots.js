/*
  Interactive ensemble slew-path and star-visit plot for an EXOSIMS summary page.

  Produces one Plotly geographic scatter plot (equirectangular projection) in
  slewPlotDiv, showing:
    - A dashed equator reference line
    - One line trace per slew pair (line width proportional to relative slew count)
    - Diamond waypoint markers near the destination of the most prominent slews,
      with hover text (Plotly cannot attach hover text to line traces on maps,
      so a nearby marker is used as a workaround)
    - Circle markers at each star (color = mean visit count across ensemble,
      size = 3 if unvisited, 8 otherwise)

  Data sources (in the directory named by ens_path_root):
    path-visits.csv  -- one row per star: name, lon, lat, visit_mean, visit_std
    path-slews.csv   -- one row per slew pair: source (index), dest (index),
                        slews (count), label

  The two CSVs are loaded in parallel via Promise.all. If either file is
  missing or empty, a red error message is written into slewPlotDiv.

  Configuration is read from data-* attributes on the slewPlotDiv element:
    data-path-root  path to the directory containing path-visits.csv and path-slews.csv
    data-path-mode  'ensemble' (default) or 'single', controls plot title and colorbar label
  No global variables or inline <script> blocks are needed.

  turmon sep 2018, revised may 2026
*/


// list of alternate names for the nearest stars, used only as a display aid
var starAlternateNames2 = {
    'HIP 5336':   ['mu Cas'],
    // automatically found, see "star-alt-names-merged.js"
    "HIP 677":    ["HD 358", "LTT 10039", "alpha And"],
    "HIP 1475":   ["HD 1326", "Gl 15A", "GX Andromedae"],
    "HIP 2021":   ["HD 2151", "Gl 19", "beta Hyi"],
    "HIP 3821":   ["HD 4614", "Gl 34A", "eta Cassiopei A"],
    "HIP 3829":   ["GJ 35", "van Maanen's Star"],
    "HIP 5643":   ["GJ 54.1", "YZ Ceti"],
    "HIP 8102":   ["HD 10700", "Gl 71", "tau Ceti"],
    "HIP 9884":   ["HD 12929", "Gl 84.3", "alf Arietis "],
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
 * Map-geometry utilities for placing hover markers partway along great-circle
 * slew arcs. Plotly cannot attach hover text to line traces on geo plots, so
 * a nearby marker is used as a workaround.
 */
function toDegrees(angle) { return angle * (180.0 / Math.PI); }
function toRadians(angle) { return angle * (Math.PI / 180.0); }

function lonlat2xyz(lon, lat) {
    return {
        x: Math.cos(toRadians(lat)) * Math.cos(toRadians(lon)),
        y: Math.cos(toRadians(lat)) * Math.sin(toRadians(lon)),
        z: Math.sin(toRadians(lat))
    };
}
function xyz2lonlat(pt) {
    var r = Math.sqrt(pt.x*pt.x + pt.y*pt.y + pt.z*pt.z);
    var lon = (toDegrees(Math.atan2(pt.y, pt.x)) + 360) % 360;
    var lat = toDegrees(Math.asin(pt.z / r));
    return {lon, lat};
}
function xyz_norm(pt) {
    return Math.sqrt(pt.x*pt.x + pt.y*pt.y + pt.z*pt.z);
}

// Return one point near the destination end of the great-circle arc from
// (lon1,lat1) to (lon2,lat2), plus r (arc chord length, 0–2).
// The constant 0.15 is empirical: sqrt(r) scaling makes it work for both
// short and long arcs.
function gcp1(lon1, lat1, lon2, lat2) {
    var p1 = lonlat2xyz(lon1, lat1);
    var p2 = lonlat2xyz(lon2, lat2);
    var delta = {x: p2.x-p1.x, y: p2.y-p1.y, z: p2.z-p1.z};
    var r = xyz_norm(delta);
    var d = 0.15 / Math.sqrt(r);
    var p_mid = {x: p2.x-d*delta.x, y: p2.y-d*delta.y, z: p2.z-d*delta.z};
    var {lon, lat} = xyz2lonlat(p_mid);
    return {r, lon, lat};
}


// Load a CSV file and parse it into an array of row objects.
function loadCSV(url) {
    return fetch(url)
        .then(function(response) {
            if (!response.ok)
                throw new Error('Could not load ' + url + ' (HTTP ' + response.status + ')');
            return response.text();
        })
        .then(function(text) {
            var lines = text.trim().split('\n');
            var headers = lines[0].split(',');
            return lines.slice(1).map(function(line) {
                var values = line.split(',');
                var row = {};
                headers.forEach(function(h, i) { row[h] = values[i]; });
                return row;
            });
        });
}

// Write a visible error message into the plot div.
function showError(div, message) {
    if (div) {
        div.innerHTML =
            '<p style="color:red; padding:20px; font-size:1.1em;">' +
            '<strong>Could not load plot data:</strong> ' + message + '</p>';
    }
}

// Return a parenthesized string of alternate star names, or '' if none known.
// Strips trailing A/B/C from the HIP name before lookup.
function formatStarAlternates(altNames, star) {
    var star_alone = star.replace(/ [A-C]$/, '');
    if (altNames.hasOwnProperty(star_alone)) {
        return ' (' + altNames[star_alone].join(', ') + ')';
    }
    return '';
}

// Build all Plotly traces from the CSV rows and render the plot.
function plotFromSlews(target, visitRows, slewRows, pathMode) {
    var slew_max  = slewRows.reduce(function(a, b) { return Math.max(a, parseFloat(b.slews)); }, 0.0);
    var thick_max = 16.0; // maximum line width in pixels

    var waypoints   = [];
    var slew_traces = [];

    // primitive equator (map objects have no real tick legends)
    slew_traces.push({
        lon: [0, 90, 180, 270, 360],
        lat: [0,  0,   0,   0,   0],
        type: 'scattergeo',
        hoverinfo: 'none',
        mode: 'lines',
        line: {color: 'rgba(0, 0, 0, 1)', dash: 'dashdot', width: 0.5}
    });

    // Accumulate slew line traces and waypoint markers
    slewRows.forEach(function(row) {
        var slews = parseFloat(row.slews);
        var label = row.label || ('Slews: ' + slews);
        if (slews < 1.0) return;

        var res = gcp1(visitRows[row.source].lon, visitRows[row.source].lat,
                       visitRows[row.dest].lon,   visitRows[row.dest].lat);
        // show a waypoint marker only for longer, prominent slews
        if (res.r > 0.15 && slews/slew_max > 0.08 && label) {
            res.text = [label,
                        'From: ' + visitRows[row.source].name,
                        'To: '   + visitRows[row.dest].name].join('<br>');
            waypoints.push(res);
        }
        slew_traces.push({
            lon: [visitRows[row.source].lon, visitRows[row.dest].lon],
            lat: [visitRows[row.source].lat, visitRows[row.dest].lat],
            type: 'scattergeo',
            // Plotly shows line tooltips only at endpoints, which conflicts with
            // the star markers, so hover is suppressed here.
            hoverinfo: 'none',
            mode: 'lines',
            opacity: 0.4,
            line: {width: (slews/slew_max) * thick_max}
        });
    });

    // Waypoint markers (diamonds near slew destinations)
    var x_way = [], y_way = [], text_way = [];
    waypoints.forEach(function(wp) {
        x_way.push(wp.lon);
        y_way.push(wp.lat);
        text_way.push(wp.text);
    });
    if (x_way.length > 0) {
        slew_traces.push({
            type: 'scattergeo',
            lon: x_way, lat: y_way,
            text: text_way,
            mode: 'markers',
            hoverinfo: 'text',
            marker: {color: 'rgba(30,30,30,0.4)', size: 7, symbol: 'diamond-open'}
        });
    }

    // Star markers (color = visit count)
    var x = [], y = [], z = [], size = [], text = [];
    visitRows.forEach(function(row) {
        var visits = parseFloat(row.visit_mean);
        var lon    = parseFloat(row.lon);
        var lat    = parseFloat(row.lat);
        x.push(lon);
        y.push(lat);
        z.push(visits === 0 ? NaN : visits);
        size.push(visits === 0 ? 3 : 8);
        text.push(row.name + formatStarAlternates(starAlternateNames2, row.name) + '<br>' +
                  visits.toFixed(3) + ' visits<br>' +
                  lon.toFixed(1) + '&deg;, ' + lat.toFixed(1) + '&deg;');
    });
    slew_traces.push({
        type: 'scattergeo',
        lon: x, lat: y,
        text: text,
        mode: 'markers',
        hoverinfo: 'text',
        marker: {
            colorscale: 'Portland',
            color: z,
            size: size,
            opacity: 0.7,
            colorbar: {
                tickfont: {family: 'Arial', size: 16, color: 'black'},
                title: {
                    text: '<br>' + (pathMode === 'single' ? 'Number of Characterizations'
                                                          : 'Mean Cumulative Characterizations'),
                    side: 'right',
                    font: {size: 16, family: 'Arial Bold'}
                }
            }
        }
    });

    insertPlotly(target, slew_traces, pathMode);
}

// Package traces and layout and render into the target div.
function insertPlotly(target, slew_traces, pathMode) {
    var plot_title = (pathMode === 'single')
        ? ['All Attempted Characterizations and Slew Paths',
           'Point Shading: Number of Visits; Slew Path Shading: Arbitrary',
           'Diamond Indicators Show Slew Information near Slew Destination'].join('<br>')
        : ['Mean Cumulative Characterizations and Slew Paths Over Ensemble',
           'Point Shading: Mean Visits across Ensemble; Slew Path Shading: Arbitrary',
           'Diamond Indicators Show Slew Information near Slew Destination'].join('<br>');

    var layout = {
        title: {text: plot_title, font: {family: 'Arial Black', size: 18, color: 'black'}},
        showlegend: false,
        geo: {
            showland: false, showcountries: false, showlakes: false, showocean: false,
            projection: {type: 'equirectangular'},
            coastlinewidth: 0,
            // axis titles do not work on geo objects in Plotly
            lataxis: {range: [-80,  80], showgrid: true},
            lonaxis: {range: [  0, 360], showgrid: true}
        },
        hovermode: 'closest'
    };
    Plotly.newPlot(target, slew_traces, layout);
}


// Entry point: read configuration from data-* attributes on the div, then load CSVs in parallel.
// DOMContentLoaded ensures the div exists whether this file is in <head> or <body>.
document.addEventListener('DOMContentLoaded', function() {
    var div = document.getElementById('slewPlotDiv');
    if (!div) return; // page has no path widget

    var pathRoot = div.dataset.pathRoot || '../path-ens';
    var pathMode = div.dataset.pathMode || 'ensemble';
    var visits_url = pathRoot + '/path-visits.csv';
    var slews_url  = pathRoot + '/path-slews.csv';

    Promise.all([loadCSV(visits_url), loadCSV(slews_url)])
        .then(function(results) {
            var visitRows = results[0];
            var slewRows  = results[1];
            if (!visitRows.length || !slewRows.length) {
                showError(div, 'One or more CSV files are empty.');
                return;
            }
            plotFromSlews('slewPlotDiv', visitRows, slewRows, pathMode);
        })
        .catch(function(err) { showError(div, err.message); });
});
