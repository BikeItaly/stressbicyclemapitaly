mapboxgl.accessToken = MAPBOX_TOKEN
var map = new mapboxgl.Map({
    container: 'map',
    style: MAPBOX_STYLE,
    zoom: 5,
    center: [12.79, 42.98],
    hash: true
});

map.addControl(new mapboxgl.NavigationControl());
map.addControl(
    new MapboxGeocoder({
        accessToken: mapboxgl.accessToken,
        mapboxgl: mapboxgl
    })
);
map.on('load', function () {
    map.addSource(STRESSCYCLELEVES_SOURCE, {
        type: 'vector',
        url: STRESSCYCLELEVES_TILESET,
        minzoom: 15,
        maxzoom: 22
    });
    map.addLayer({
        'id': STRESSCYCLELEVES_LAYER,
        'type': 'line',
        'source': STRESSCYCLELEVES_SOURCE,
        'source-layer': STRESSCYCLELEVES_LAYER,
        'layout': {
            'line-join': 'round',
            'line-cap': 'round'
        },
        'paint': {
            'line-width': 3,
            'line-color': ['match', ['string', ['get', 'level']], ...LEVELSANDCOLORS, '#AAAAAA'],
            'line-opacity': 1,
        },
    });
});

map.on('load', function () {
    map.addSource(STRESSCYCLELEVES_SOURCEmin, {
        'type': 'vector',
        'tiles': [
            STRESSCYCLELEVES_TILESETmin
        ],
        'minzoom': 5,
        'maxzoom': 14
    });
    map.addLayer(
        {
            'id': STRESSCYCLELEVES_LAYERmin,
            'type': 'line',
            'source': STRESSCYCLELEVES_SOURCEmin,
            'source-layer': STRESSCYCLELEVES_LAYERmin,
            'layout': {
                'line-cap': 'round',
                'line-join': 'round'
            },
            'paint': {
                'line-color': ['match', ['string', ['get', 'level']], ...LEVELSANDCOLORS, '#AAAAAA'],
                'line-width': 2
            }
        }
    );
});



function getLevelsAndColors() {
    ret = [];
    for (var k in levels) {
        ret.push(k);
        ret.push(levels[k].visible ? levels[k].color : 'rgba(0,0,0,0)');
    }
    return ret;
}

var filterGroup = document.getElementById('filter-group');

["level_1", "level_2", "level_3", "level_4"].forEach(function (layerID) {
    var input = document.createElement('input');
    input.type = 'checkbox';
    input.id = layerID;
    input.checked = true;
    input.style.backgroundColor = levels[layerID].color;
    filterGroup.appendChild(input);

    var label = document.createElement('label');
    label.setAttribute('for', layerID);
    label.textContent = levels[layerID].label;
    label.style.backgroundColor = levels[layerID].color;
    filterGroup.appendChild(label);


    // When the checkbox changes, update the visibility of the layer.
    input.addEventListener('change', function (e) {
        levels[layerID].visible = e.target.checked;

        var LevelsAndColors = getLevelsAndColors();
        map.setPaintProperty(STRESSCYCLELEVES_LAYER, 'line-color', ['match', ['string', ['get', 'level']], ...LevelsAndColors, '#AAAAAA']);
        map.setPaintProperty(STRESSCYCLELEVES_LAYERmin, 'line-color', ['match', ['string', ['get', 'level']], ...LevelsAndColors, '#AAAAAA']);
    });
});
