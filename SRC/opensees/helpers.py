
def node_average(model, response, ndm=2, keys=None):
    """
    Return a dictionary of stress components per node.
    For 2D: sxx, syy, sxy
    For 3D: sxx, syy, szz, sxy, syz, sxz
    Summed values from each element's node contribution are averaged at the end.

    Args:
        model: an OpenSees model-like object with the following interface:
               - getEleTags()
               - getEleClassTags(eleTag)
               - getNodeTags()
               - eleNodes(eleTag)
               - eleResponse(eleTag, str)
        response (str): a label if needed for e.g. 'stress', 'stressAtNodes'
        ndm (int): number of dimensions (2 or 3)

    Returns:
        dict[nodeTag] = {
            "xx": float, 
            "yy": float,
            (optionally in 2D) "xy": float,
            (optionally in 3D) "zz", "yz", "xz",
            ... possibly derived quantities ...
        }
    """

    # Helper: pick which elements are 2D vs 3D, if needed
    # For example, you might have separate function(s) or logic:
    #    ele_tags = _stress_2d_ele_tags_only(model, model.getEleTags())
    # or for 3D:
    #    ele_tags = _stress_3d_ele_tags_only(model, model.getEleTags())
    #
    # Below, we simply re-use the user's helper for 2D, but you'll have to define
    # your own for 3D or unify them if your model can contain both:
    ele_tags_all = model.getEleTags()
    if ndm == 2:
        ele_tags = _stress_2d_ele_tags_only(model, ele_tags_all)
    else:
        # Placeholder; you must implement a 3D filter or skip if you already know
        # the elements are 3D:
        ele_tags = ele_tags_all

    # For a typical 2D element with 3 stress components: sxx, syy, sxy
    # For a typical 3D element with 6 stress components: sxx, syy, szz, sxy, syz, sxz
    # You may need to adapt this depending on your element type(s).
    if ndm == 2:
        keys = "sxx", "syy", "sxy"
    else:
        keys = "sxx", "syy", "szz", "sxy", "syz", "sxz"


    node_tags = model.getNodeTags()

    # Initialize a dictionary to hold stress sums and a count of how many
    # times each node is encountered (to do averaging).
    sig_dict = {node: {"_count": 0, **{key: 0.0 for key in keys}} for node in node_tags}

    nrc = len(keys)

    # We assume each element’s node-level stress output is available:
    # e.g. model.eleResponse(elem_tag, 'stressAtNodes')
    # and yields a flat array (with length = n_nodes_in_element * nrc).
    # Then we reshape it to (nen, nrc).
    for elem in ele_tags:
        ele_node_tags = model.eleNodes(elem)
        nen = len(ele_node_tags)

        sig_per_node = np.reshape(model.eleResponse(elem, response), (nen, nrc))

        # Accumulate into sig_dict
        for i, node in enumerate(ele_node_tags):
            for j, key in enumerate(keys):
                sig_dict[node][key] += sig_per_node[i][j]

            sig_dict[node]["_count"] += 1

    # Now do averaging (if a node belongs to multiple elements).
    for node in node_tags:
        c = sig_dict[node]["_count"] or 1.0

        for key in keys:
            sig_dict[node][key] /= c


    # Clean up
    for node in node_tags:
        sig_dict[node].pop("_count", None)

    return sig_dict


def find_node(model, x=None, y=None, tol=1e-8):
    """
    Return the tag of the *first* node that matches the specified coordinates
    within a given tolerance. If no match is found, returns None.

    Parameters
    ----------
    model : opensees.openseespy.Model
        The model object containing nodes.
    x : float, optional
        X-coordinate to match (if None, no X match is required).
    y : float, optional
        Y-coordinate to match (if None, no Y match is required).
    tol : float
        Matching tolerance for coordinates.

    Returns
    -------
    int or None
        The tag of the first matching node, or None if not found.
    """
    # Get the list of all node tags in the model
    node_tags = model.getNodeTags()

    for tag in node_tags:
        coords = model.nodeCoord(tag)  # [xCoord, yCoord, ...]
        # Check X coordinate if requested
        if x is not None and abs(coords[0] - x) > tol:
            continue
        # Check Y coordinate if requested
        if y is not None and abs(coords[1] - y) > tol:
            continue

        # If we've passed all checks, this node matches
        return tag

    # If no node passes the test, return None
    return None


def find_nodes(model, x=None, y=None, z=None, tol=1e-8):
    """
    Return an iterator of node tags that match the specified coordinates
    within a given tolerance.

    Parameters
    ----------
    model : opensees.openseespy.Model
        The model object containing nodes.
    x : float, optional
        X-coordinate to match (if None, no X match is required).
    y : float, optional
        Y-coordinate to match (if None, no Y match is required).
    tol : float
        Matching tolerance for coordinates.

    Yields
    ------
    int
        The tag of each matching node.
    """
    node_tags = model.getNodeTags()
    for tag in node_tags:
        coords = model.nodeCoord(tag)  # [xCoord, yCoord, ...]
        # Check X coordinate if requested
        if x is not None and abs(coords[0] - x) > tol:
            continue
        # Check Y coordinate if requested
        if y is not None and abs(coords[1] - y) > tol:
            continue
        # Check Z coordinate if requested
        if z is not None and abs(coords[2] - z) > tol:
            continue
        # If we've passed all checks, this node matches
        yield tag

