mapboxgl.accessToken = MAPBOX_TOKEN
var map = new mapboxgl.Map({
    container: 'map',
    style: MAPBOX_STYLE,
    zoom: 15,
    center: [11.1234, 46.0686],
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
        minzoom: 12,
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