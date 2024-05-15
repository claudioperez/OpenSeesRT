from .section import (
      patch,
      layer,
      FiberSection,
      ConfinedPolygon,
      ConfiningPolygon,
      PolygonRing
)

SectionGeometry = FiberSection

def _WideFlange(material, tag=None, ndm=3, **sec_data)->SectionGeometry:

    d  = sec_data['d']
    bf = sec_data['bf']
    tf = sec_data['tf']
    tw = sec_data['tw']

    nft = sec_data['nft']
    nwl = sec_data['nwl']
    nfl = sec_data.get('nfl', 1, ) #sec_data['nft'])
    nwt = sec_data.get('nwt', 1, ) #sec_data['nwl'])

    int_typ = sec_data.get("IntTyp", None)
    flg_opt = sec_data.get('FlgOpt', True)
    web_opt = sec_data.get('WebOpt', False)

    yoff = ( d - tf) / 2
    zoff = (bf + tw) / 4

    dw = d - 2 * tf
    bi = bf - tw

    if ndm == 2:
        GJ =  1.0

    else:
        J = kwds.get("J")

        # return SectionGeometry(name=tag, GJ=GJ, areas=[
        #     patch.rect(vertices=[[-],[]], material=material)
        #     # patch.Fiber([x,y], A, material) for x,y,A in zip(xfib, yfib, wfib)
        # ])

    if flg_opt:
        return SectionGeometry(name=tag, GJ=GJ, areas=[
            patch.rect(corners=[[-bf/2, yoff-tf/2],[bf/2,  yoff+tf/2]], material=material, divs=(nfl, nft), rule=int_typ),
            patch.rect(corners=[[-tw/2,-yoff+tf/2],[tw/2,  yoff-tf/2]], material=material, divs=(nwt, nwl), rule=int_typ),
            patch.rect(corners=[[-bf/2,-yoff-tf/2],[bf/2, -yoff+tf/2]], material=material, divs=(nfl, nft), rule=int_typ),
        ])

    else:
        return SectionGeometry(name=tag, GJ=GJ, areas=[
            patch.rect(corners=[[-zoff-bi/4, yoff-tf/2],[-zoff+bi/4,  yoff+tf/2]], material=material, divs=(nfl, nft), rule=int_typ),
            patch.rect(corners=[[ zoff-bi/4, yoff-tf/2],[ zoff+bi/4,  yoff+tf/2]], material=material, divs=(nfl, nft), rule=int_typ),
            patch.rect(corners=[[     -tw/2,-yoff-tf/2],[      tw/2,  yoff+tf/2]], material=material, divs=(nwt, nwl), rule=int_typ),
            patch.rect(corners=[[-zoff-bi/4,-yoff-tf/2],[-zoff+bi/4, -yoff+tf/2]], material=material, divs=(nfl, nft), rule=int_typ),
            patch.rect(corners=[[ zoff-bi/4,-yoff-tf/2],[ zoff+bi/4, -yoff+tf/2]], material=material, divs=(nfl, nft), rule=int_typ),
        ])



def from_aisc(type, identifier, material = None, tag: int = None, mesh:dict=None, units=None, ndm=3, **kwds):
    if mesh is None:
        mesh = {}
    if units is None:
        import opensees.units.english as units

    aisc_data = load_aisc(identifier, units=units)


    if identifier[0] == "W":
        geom = _WideFlange(**aisc_data, **mesh, tag=tag, ndm=ndm, material=material)


    if type == "Fiber":
        return geom
    if type == "Resultant":
        pass
    if type == "Elastic":
        pass

def load_aisc(SectionName, props="", units=None)->dict:
    """Load cross section properties from AISC database.

    props:
        A list of AISC properties, or one of the following:
        - 'simple': `A`, `Ix`, `Zx`

    """

    from .aisc_imperial import imperial
    SectData = imperial[SectionName.upper()]
    if props == "simple":
        props = ""
        return

    elif props:
        props = props.replace(" ", "").split(",")
        sectData = {k: v for k, v in SectData.items() if k in props}
        if "I" in props:
            sectData.update({"I": SectData["Ix"]})
        return sectData
    
    for k,v in list(SectData.items()):
        try:
            SectData[k] = float(v)
        except:
            pass

    return SectData
