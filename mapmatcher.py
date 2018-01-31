# Code adapted an debugged for this research project purposes. September 2017

# Author:   Nimol Vamoeurn 2017
#           nimolva@gmail.com
# Purpose:  This script map matches point sequences to an underlying street
#           network using the Hidden Markov Model and the Viterbi algorithm
#           with probability parameterized based on spatial and network
#           distances.
# Credit:   Functions in the original code from the work of Nimol Vamoeurn
#           were used and adapted for the map matching processes.
#           The code was originally adapted from the work of
#           Simon Scheider (2017) available at
#           https://github.com/simonscheider/mapmatching



from shapely.geometry import shape, Point, mapping
from shapely.ops import linemerge
import fiona
import networkx as nx
import csv
from rtree import index
from math import exp, sqrt
from datetime import datetime
import itertools
from collections import Counter, defaultdict

startTime = datetime.now()


# Main map matching mathod based on the Viterbi algorithm for the Hidden Markov
# Model
# points = sequence of tracking data ordered by time
# segmentInfo = refer to getSegmentInfo method
# graph = graph of the street network (refer to getNetworkInfo method)
# tree = R-Tree spatial index (refer to buildRTree method)
# segmentShapes = binary geom of street segments
# decayConstantNet = network distance (meters) related to how far track point
# are separated from one another
# decayConstantEu = Euclidean distance (meters) related to positional
# accuaracy of track points (accuracy of GPS receiver)
# maxDist = Euclidean distance as limit farther than which segment
# candidates are ignored
def mapMatch(points, segmentInfo, graph, tree, segmentShapes,
             decayconstantNet=30, decayConstantEu=10, maxDist=50):

    # Ensure that distance parameters are floats
    decayconstantNet = float(decayconstantNet)
    decayConstantEu = float(decayConstantEu)
    maxDist= float(maxDist)

    # Probability distribution (based on Viterbi algorithm) of street
    # segments for a given GPS track point, with the most probable predecessor
    #  segment taken into account
    V = [{}]

    endpoints = segmentInfo[0]  # endpoint coordinates of all street segments
    lengths = segmentInfo[1]    # length of all street segments
    pathnodes = []              # set of path nodes to prevent loops

    # Get segment candidates for the first point of the track
    sc = getSegmentCandidates(points[0], tree, segmentShapes, decayConstantEu,maxDist)

    for s in sc:
        V[0][s] = {"prob": sc[s], "prev": None, "path": [], "pathnodes": []}

    #print V
    # Run Viterbi algorithm when t > 0
    for t in range(1, len(points)):
        V.append({})
        # Store previous segment candidates
        lastsc = sc
        # Get segment candidates and their a priori probability
        sc = getSegmentCandidates(points[t], tree, segmentShapes,
                                  decayConstantEu, maxDist)
        for s in sc:
            max_tr_prob = 0     # init maximum transition probability
            prev_ss = None
            path = []
            for prev_s in lastsc:
                # determine the highest transition probability from previous
                # candidates to s and get the corresponding network path
                pathnodes = V[t-1][prev_s]["pathnodes"][-10:]
                n = getNetworkTransP(prev_s, s, graph, endpoints,
                                     lengths, pathnodes, decayconstantNet)
                np = n[0]   # network transition probability
                tr_prob = V[t-1][prev_s]["prob"]*np
                # Select the most probable predecessor candidate and the
                # path to it
                if tr_prob > max_tr_prob:
                    max_tr_prob = tr_prob
                    prev_ss = prev_s
                    path = n[1]
                    if n[2] is not None:
                        pathnodes.append(n[2])
            # Final probability of a candidate is the product of a-priori
            # and network transitional probability
            max_prob = sc[s] * max_tr_prob
            V[t][s] = {"prob": max_prob, "prev": prev_ss, "path": path,
                       "pathnodes": pathnodes}

        maxv = max(value["prob"] for value in V[t].values())
        maxv = (1 if maxv == 0 else maxv)
        for s in V[t].keys():
            V[t][s]["prob"] = V[t][s]["prob"]/maxv

    # opt is the result: a list of matched segments [s1, s2, s3,...] in
    # the exact order of the track points: [p1, p2, p3,...]
    opt = []

    # Get the highest probability at the end of the track
    max_prob = max(value["prob"] for value in V[-1].values())
    previous = None
    #if max_prob == 0:
        #print(" probabilities fall to zero (network distances in data are " \
             # "too large, try increasing network decay parameter)")

    # Get the most probable ending state and its backtrack
    for st, data in V[-1].items():
        if data["prob"] == max_prob:
            opt.append(st)
            previous = st
            break

    # Follow the backtrack till the first observation to fish out most
    # probable states and corresponding paths
    for t in range(len(V) - 2, -1, -1):
        # Get the subpath between last and most probable previous segment and
        # add it to the resulting path
        path = V[t + 1][previous]["path"]
        opt[0:0] = (path if path != None else [])
        # Insert the previous segment
        opt.insert(0, V[t + 1][previous]["prev"])
        previous = V[t + 1][previous]["prev"]

    # Clean the path (remove double segments and crossings, etc.)
    #opt = cleanPath(opt, endpoints)
    return opt


# Clean the path by removing redundant or unnecessary segments from the path
def cleanPath(opt, endpoints):
    last =()
    lastlast =()
    optout = []
    for s in opt:
        if s != last:
            match = False
            if last != () and lastlast != ():
                lastep = endpoints[last]
                lastlastep = endpoints[lastlast]
                sep = endpoints[s]
                for j in lastlastep:
                    if lastep[0]== j:
                        for k in sep:
                            if lastep[1] == k:
                                match = True
                    elif lastep[1]== j:
                        for k in sep:
                            if lastep[0] == k:
                                match = True
            elif last != ():
                sep = endpoints[s]
                lastep = endpoints[last]
                for k in sep:
                    if lastep[1] == k or lastep[0] == k:
                        match = True
            if match:
                optout.append(last)
            if s == opt[-1]:
                optout.append(s)
        lastlast = last
        last = s
    return optout


# Export all paths from different tracking sessions into a single shapefile
# all_opts: a dictionary where key is tracking session id and value is a
# list of all street segments traversed by the tracking session
# shp_out: output shapefile to be exported
# segments: the shapefile representing the street network
def exportMultiPaths(all_opts, shp_out, segments):
    with fiona.open(segments, 'r') as source:
        source_crs = source.crs
        new_schema = {'geometry': 'LineString',
                      'properties': {'ses_id': 'str:35'}}
        with fiona.open(shp_out, 'w', driver='ESRI Shapefile',
                        crs=source_crs, schema=new_schema) as track:
            for key, opt in all_opts.items():
                # Discard empty paths corresponding to tracking points that
                # do not produce a meaningful movement
                if len(opt) > 1:
                    path = []
                    prop_t = {'ses_id': key}
                    for i in opt:
                        geom = shape(source[i]['geometry'])
                        path.append(geom)
                    track.write({'properties': prop_t, 'geometry': mapping(
                        linemerge(path))})
                else:
                    pass


# Export tracking sessions that do not produce any meaningful movement into
# a CSV file
def exportEmptyPaths(all_opts, csv_out):
    with open(csv_out, 'wb') as csv_file:
        writer = csv.writer(csv_file)
        for key, opt in all_opts.items():
            if len(opt) <= 1:
                writer.writerow([key])


# Export all path into a CSV file where first column shows the street
# segment and the second column shows its frequency of being used by all paths
def exportFreqCSV(all_opts, csv_out):
    with open(csv_out, 'wb') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['Segment_ID', 'Frequency'])
        opts_lst = []
        for key, opt in all_opts.items():
            if len(opt) > 1:
                opts_lst.append(opt)
        merged_opts = list(itertools.chain.from_iterable(opts_lst))
        segment_freq = dict(Counter(merged_opts))
        for key, value in segment_freq.items():
            writer.writerow([key, value])


# Export a CSV file whose first column represents street segment ID and
# each subsequent column representing the session ID of the track that pass
# through that segment
def exportSegmentPaths(all_opts, csv_out):
    with open(csv_out, 'wb') as csv_file:
        writer = csv.writer(csv_file)
        segment_paths = defaultdict(set)
        for key, values in all_opts.items():
            for value in values:
                segment_paths[value].add(key)

        for segment, paths in segment_paths.items():
            writer.writerow([segment,paths])


# Export a single tracking session as a single shapefile containing only one
# line feature representing the path of the track
# trackname: name of the exported shapefile
def exportPath(opt, trackname, segments):
    with fiona.open(segments, 'r') as source:
        source_crs = source.crs
        source_schema = source.schema
        with fiona.open(trackname, 'w', driver='ESRI Shapefile',
                        crs=source_crs, schema=source_schema) as track:
            path = []
            propt = source[opt[0]]['properties']
            for i in opt:
                geom = shape(source[i]['geometry'])
                path.append(geom)
            track.write({'properties': propt, 'geometry': mapping(
                linemerge(path))})


# Export a single tracking session as a single shapefile containing multiple
# line features each of which represents the street segment of the track's path
# trackname: name of the exported shapefile
def exportPathSegments(opt, trackname, segments):
    with fiona.open(segments, 'r') as source:
        source_crs = source.crs
        source_schema = source.schema
        with fiona.open(trackname, 'w', driver='ESRI Shapefile',
                        crs=source_crs, schema=source_schema) as track:
            for i in opt:
                propt = source[i]['properties']
                geom = shape(source[i]['geometry'])
                track.write({'properties': propt, 'geometry': mapping(
                    geom)})


# Return closest segment candidates with a-priori probabilities based on
# maximal spatial distance of segments from point
def getSegmentCandidates(point, tree, segmentShapes, decayConstantEu,
                         maxdist=50):
    candidates = {}
    make_point = Point(point)
    point_buffer = make_point.buffer(maxdist)
    for f_id in list(tree.intersection(point_buffer.bounds)):
        line = segmentShapes.get(f_id)
        dist = make_point.distance(line)
        if dist <= maxdist:
            candidates[f_id] = getPDProbability(dist, decayConstantEu)
    return candidates


# Return transition probability of going from segment s1 to s2
def getNetworkTransP(s1, s2, graph, endpoints, segmentlengths, pathnodes,
                     decayconstantNet):
    subpath = []
    s1_point = None
    s2_point = None

    if s1 == s2:
        dist = 0
    else:
        # Obtain edges (tuples of endpoints) for segment identifiers
        s1_edge = endpoints[s1]
        s2_edge = endpoints[s2]

        # Determines segment endpoints of the two segments that are
        # closest to each other
        minpair = [0, 0, 100000]
        for i in range(0, 2):
            for j in range(0, 2):
                d = round(pointdistance(s1_edge[i], s2_edge[j]), 2)
                if d < minpair[2]:
                    minpair = [i, j, d]
        s1_point = s1_edge[minpair[0]]
        s2_point = s2_edge[minpair[1]]

        if s1_point == s2_point:
            # If segments are touching, use a small network distance
            dist = 5
        else:
            try:
                # Compute the shortest path based on segment length on
                # street network graph where segment endpoints are nodes and
                #  segments are (undirected) edges
                if graph.has_node(s1_point) and graph.has_node(s2_point):
                    dist = nx.astar_path_length(graph, s1_point, s2_point,
                                                weight='length',
                                                heuristic=pointdistance)
                    path = nx.astar_path(graph, s1_point, s2_point,
                                         weight='length',
                                         heuristic=pointdistance)
                    path_edges = zip(path, path[1:])
                    subpath = []
                    for e in path_edges:
                        oid = graph.edge[e[0]][e[1]]["OBJECTID"]
                        subpath.append(oid)
                else:
                    dist = 3*decayconstantNet
            except nx.NetworkXNoPath:
                dist = 3*decayconstantNet
    return (getNDProbability(dist, decayconstantNet), subpath, s2_point)


# Return the probability that a point is on a segment
# dist: Euclidean distance between the point and the segment (in meter)
# decayconstant: Euclidean distance in meter farther than which probability
# falls to 0.34
def getPDProbability(dist, decayconstant=10):
    decayconstant = float(decayconstant)
    dist= float(dist)
    try:
        p = 1 if dist == 0 else round(1/exp(dist/decayconstant), 4)
    except OverflowError:
        p = round(1/float('inf'), 2)
    return p


# Return the probability that a segment is the successor of another on a track
# dist: Euclidean distance between two segments (in meter)
# decayconstant: Network distance in meter farther than which probability
# falls to 0.34
def getNDProbability(dist, decayconstant=30):
    decayconstant = float(decayconstant)
    dist = float(dist)
    try:
        p = 1 if dist == 0 else round(1/exp(dist/decayconstant), 2)
    except OverflowError:
        p = round(1/float('inf'), 2)
    return p


# Build an RTREE spatial index of the street network
# shapefile: a shapefile containing street segments (must be planarized)
def buildRTree(shapefile):
    idx = index.Index()
    with fiona.open(shapefile, 'r') as streets:
        for st in streets:
            st_id = int(st['properties']['OBJECTID'])
            st_geom = shape(st['geometry'])
            idx.insert(st_id, st_geom.bounds)
    return idx


# Build a network graph of the street
# shapefile: a shapefile containing street segments (must be planarized)
# segmentlengths: obtained from getSegmentInfo
def getNetworkGraph(shapefile, segmentlengths):
    g = nx.read_shp(shapefile)
    sg = g.to_undirected()
    for n0, n1 in sg.edges_iter():
        oid = int(sg[n0][n1]["OBJECTID"])
        sg.edge[n0][n1]['length'] = segmentlengths[oid]
    return sg


# Returns a dictionary containing OBJECTID key and binary geom value
# shapefile: a shapefile containing street segments (must be planarized)
def getSegmentShapes(shapefile):
    shapes = {}
    with fiona.open(shapefile, 'r') as streets:
        for st in streets:
            st_id = int(st['properties']['OBJECTID'])
            st_geom = shape(st['geometry'])
            shapes[st_id] = st_geom
    return shapes


# Returns endpoints coordinates of all segments and their length
# shapefile: a shapefile containing street segments (must be planarized)
def getSegmentInfo(shapefile):
    with fiona.open(shapefile, 'r') as streets:
        endpoints = {}
        segmentlength = {}
        for segment in streets:
            seg_id = int(segment['properties']['OBJECTID'])
            seg_coord = segment['geometry']['coordinates']
            seg_length = shape(segment['geometry']).length
            endpoints[seg_id] = (seg_coord[0], seg_coord[-1])
            segmentlength[seg_id] = seg_length
        return endpoints, segmentlength


# Return a list containing coordinates of all points stored in a shapefile
# points_shp: a shapefile containing track points
def getPointsSHP(points_shp):
    with fiona.open(points_shp, 'r') as input_track:
        trackpoints = []
        for point in input_track:
            point_coord = point['geometry']['coordinates']
            trackpoints.append(point_coord)
        return trackpoints


# Return a list containing coordinates of all points stored in a CSV file
# points_csv: a CSV file containing track points
def getPointsCSV(points_csv):
    with open(points_csv, 'r') as input_track:
        track = csv.reader(input_track, delimiter=',')
        trackpoints = []
        for point in track:
            trackpoints.append((int(point[1]), int(point[2])))
        return trackpoints


# Return distance between two points p1 and p2
# p1 or p2 must be a tuple of x and y coordinates
def pointdistance(p1, p2):
    dist = sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)
    return dist
