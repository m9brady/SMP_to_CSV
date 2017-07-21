try:
    import os
    import numpy as np
    import geojson, folium
    from folium.plugins import MarkerCluster, HeatMap
except ImportError:
    print("Unable to import required modules. Did you forget to install folium, geojson?")


def export_site_map(smpList, outHtmlFile, center=(39.038759, -108.015841), mapType='Cluster'):
    '''
    Using the folium module, create and save a browser-renderable page displaying the SMP sites with a terrain basemap.
    
    Parameters
    ----------
    smpList : 
        The list of SMP objects to be included in the map.
    outHtmlFile : 
        The absolute path to the output html file.
    center : (float, float)     (default = (39.038759, -108.015841))
        The map center point in ``(lat, lon)`` format. Default is Grand Mesa, CO, USA.
    mapType : {'Cluster' or 'Heat'}
        Can be one of either ``Cluster`` for Leaflet MarkerCluster-styled points, or ``Heat`` to generate a heatmap rather than point features.
    
    '''
    geomList = []
    for s in smpList:
        siteName = s.header['File Name']
        coords = [s.header['Longitude'], s.header['Latitude']]
        geom = geojson.Feature(siteName, geojson.Point(coords))
        geomList.append(geom)
    #featColl = geojson.FeatureCollection(geomList)
    lats = np.array([g.geometry.coordinates[1] for g in geomList])
    lons = np.array([g.geometry.coordinates[0] for g in geomList])
    fMap = folium.Map(location=center, 
                   tiles='Stamen Terrain',
                   zoom_start=12)
    #urCorner = [max(lats+90)-90, max(lons+180)-180]
    #llCorner = [min(lats+90)-90, min(lons+180)-180]
    #fMap.fit_bounds([llCorner, urCorner])
    
    if mapType.lower() == 'cluster':
        popups = [g.id for g in geomList]
        fMap.add_child(MarkerCluster(locations=list(zip(lats,lons)), popups=popups))
    elif mapType.lower() == 'heat':
        fMap.add_child(HeatMap(zip(lats, lons, np.repeat(1, len(lats))), radius=10))
    try:
        fMap.save(outHtmlFile)
    except:
        print("Unable to save file: {}".format(os.path.abspath(outHtmlFile)))
    