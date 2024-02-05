import numpy as np

def clean_line(line):
    # Some integrators do produce lines with sequences of multiple identical points
    i, k = 0, line.shape[0]-1
    while i < k:
        if np.array_equal(line[i+1], line[i]):
            line = np.delete(line, i, axis=0)
            k -= 1
        else:
            i += 1
    return(line)

def get_cumulative_angles(line):
    # Calculate vectors between consecutive points
    vec = line[1:]-line[:-1]
    # Calculate angles between consecutive vectors
    (vx1, vy1), (vx2, vy2) = vec[:-1].T, vec[1:].T
    angles = np.arctan2(vx1*vy2-vx2*vy1, vx1*vx2+vy1*vy2)
    # Cumulative sum of angles
    return np.cumsum(angles)

def has_loop(line):
    cumulative_angles = get_cumulative_angles(line)
    # Loops are detected after a + or - 2*pi rotation, that heuristic may fail in some cases
    if abs(cumulative_angles[-1]) < 2*np.pi: return 0  # There is no loop
    if cumulative_angles[-1] > +2*np.pi: return +1 # Trigonometric loop
    if cumulative_angles[-1] < -2*np.pi: return -1 # Clockwise loop

def projection_lambda(line, point, i):
    p0, p1 = line[i:i+2]
    return np.sum((point-p0)*(p1-p0))/np.sum((p1-p0)*(p1-p0))

def get_quasi_loop(line):
    # Quasi-loop is line interrupted at the projection of the seed point after a + or - 2*pi rotation
    cumulative_angles = get_cumulative_angles(line)
    # Start after approximately a + or - 2*pi rotation, that heuristic may fail in some cases
    i = np.argmax(abs(cumulative_angles) > 2*np.pi)
    # Find the projection of the seed
    previous_move = 0  # Avoid infinite back and forth if the projection is on a corner
    while True:
        l = projection_lambda(line, line[0], i)
        if l < 0:
            if previous_move == +1 or i+2 == len(line):
                break  # Projection is on a corner or we would get out of the line
            else:
                i -= 1
                previous_move = -1
        elif l > 1:
            if previous_move == -1 or i == 0:
                break  # Projection is on a corner or we would get out of the line
            else:
                i += 1
                previous_move = +1
        else:
            break  # Projection of teh seed is within the i to i+1 segment, bounds included
    if i == 0 or i+2 >= len(line):
        return None  # We should not be at an end of the line
    l = np.clip(projection_lambda(line, line[0], i), 0, 1)
    quasi_loop = line[0:i+1]
    if l > 0:
        quasi_loop = np.concatenate([quasi_loop, [(1-l)*line[i]+ l*line[i+1]]], axis=0)
    return clean_line(quasi_loop)

def close_loop(line, max_iter=50, verbose=0):
    previous_drift, stop = np.inf, False
    for k in range(max_iter):
        # print(k, previous_drift)
        quasi_loop = get_quasi_loop(line)
        if quasi_loop is None:
            return line, np.inf
        # Quasi-loop normal at the end
        vectors = quasi_loop[1:]-quasi_loop[:-1]
        normal = np.array([vectors[-1,1], -vectors[-1,0]])
        normal /= np.linalg.norm(normal)
        # Quasi-loop drift
        drift = np.sum(normal*(quasi_loop[-1]-quasi_loop[0]))
        if previous_drift == np.inf:
            initial_drift = drift
        if abs(drift) < abs(previous_drift):
            previous_drift = drift
            previous_quasi_loop = quasi_loop
        else:
            break
        # Quasi-loop length
        lengths = np.linalg.norm(vectors, axis=1)
        length = np.sum(lengths)
        # Quasi-loop center
        centers = 0.5*(quasi_loop[:-1]+quasi_loop[1:])
        center = np.sum(centers*np.expand_dims(lengths, axis=1), axis=0)/length
        # From center unit vectors
        nr, nrpoint = line-center, quasi_loop[-1]-center
        nrnorm, nrpointnorm = np.linalg.norm(nr,axis=1, keepdims=True), np.linalg.norm(nrpoint)
        nrnorm, nrpointnorm = np.where(nrnorm == 0, 1, nrnorm), 1 if nrpointnorm == 0 else 1
        nr, nrpoint = nr/nrnorm, nrpoint/nrpointnorm
        # Drift correction factor
        dcf = drift/length/np.sum(normal*nrpoint)
        # Correct drift
        padded_cumulative_lengths = np.pad(np.cumsum(np.linalg.norm(line[1:]-line[:-1], axis=1)), (1, 0))
        line = line-dcf*nr*np.expand_dims(padded_cumulative_lengths, axis=1)
    if verbose > 0:
        print(k+1, "iterations, initial drift =", initial_drift, "final_drift =", previous_drift)
    previous_quasi_loop[-1] = previous_quasi_loop[0]  # Actually close the loop after drift is minimized
    return previous_quasi_loop, initial_drift

def merge_or_close_pair(pair, max_iter=50, verbose=0):
    # Assumes this order from streamlines_from_source(), NOT checked
    fline, bline = pair 
    fline, bline = clean_line(fline), clean_line(bline)
    floop, bloop = has_loop(fline), has_loop(bline)
    if floop == 0 and bloop == 0: # No loop
        if verbose > 0:
            print("No loop")
        return np.concatenate([np.flip(bline, axis=0), fline[1:]])
    if floop != 0 and bloop == 0: # Loop on forward line only, infrequent, not tested
        if verbose > 0:
            print("Loop on forward line only")
        return close_loop(fline)[0]
    if floop == 0 and bloop != 0: # Loop on forward line only, infrequent, not tested
        if verbose > 0:
            print("Loop on backward line only")
        return np.flip(close_loop(bline)[0], axis=0)
    # Loop on both forward and backward loops, choose the one with the smallest initial drift
    if verbose > 0:
        print("Loop on both forward and backward lines")
    fline, fdrift = close_loop(fline, max_iter=max_iter, verbose=verbose)
    bline, bdrift = close_loop(bline, max_iter=max_iter, verbose=verbose)
    if fdrift == np.inf and bdrift == np.inf:
        if verbose > 0:
            print("failed to close loop")
        return fline [0:1]
    if abs(fdrift) < abs(bdrift):
        if verbose > 0:
            print("forward line has the smallest initial drift")
        return fline
    else:
        if verbose > 0:
            print("backward line has the smallest initial drift")
        return np.flip(bline, axis=0)
    
def merge_or_close_pairs(pairs, max_iter=50, verbose=0):
    plines = []
    for i, pair in enumerate(pairs):
        if len(pair) == 2:
            if verbose > 0:
                print("Seed", i, end = " ")
            plines.append(merge_or_close_pair(pair, verbose=verbose))
        elif verbose >= 0:
            print("Seed", i, "p : forward or backward line is missing,", len(pair), "lines")
    return plines

# For use with pyvista
# Extract pairs of raw streamlines associated to seeds from a pyvista streamlines object
def get_pv_pairs(streamlines, seeds):
    # get the lines and points arrays
    lines = streamlines.lines
    points = streamlines.points
    # iterate over the lines array and extract the point coordinates
    line_points = []
    i = 0
    while i < lines.shape[0]:
        l = lines[i]
        inds = lines[i+1:i+1+l]
        i += l+1
        pts = points[inds][:,:2]
        line_points.append(pts)
    # Assumes [forward, backward] order
    pairs = [[line for line in line_points if np.array_equal(line[0], point[0:2])] for point in seeds.points]
    return pairs

# For use with matplotlib
# Extract raw line pairs associated to seeds from a matplotlib streamplot object
def get_mpl_pairs(splt, seeds):
    segments = splt.lines.get_segments()
    if len(segments) == len(seeds): # Format in 3.8.0, list of lines as arrays of points
        # Remove duplicate points
        lines = []
        for line_points in segments:
            line = [point for i, point in enumerate(line_points) if i == 0 or not np.array_equal(point, line_points[i-1])]
            lines.append(np.array(line))
    else: # Format in matplotlib 3.7.2, all segments for all lines in a single list
        # Remove null segments
        non_nulls = []
        for i in range(len(segments)):
            if not np.array_equal(segments[i][0], segments[i][1]):
                non_nulls.append(segments[i])
        segments = non_nulls
        # Split into lines
        breaks = []
        for i in range(len(segments)-1):
            if not np.array_equal(segments[i+1][0], segments[i][1]):
                breaks.append(i+1)
        lines_segments = [segments[:breaks[0]]]
        for i in range(len(breaks)-1):
            lines_segments.append(segments[breaks[i]:breaks[i+1]])
        lines_segments.append(segments[breaks[-1]:])
        # Turn array of segments into array of points, line by line
        lines = []
        for line_segments in lines_segments:
            line = [line_segments[0][0]]
            for segment in line_segments:
                line.append(segment[1])
            lines.append(np.array(line))
    # Split lines into [forward, backward] pairs, assume lines goes downstream
    pairs = []
    for line, seed in zip(lines, seeds):
        # Assume that the lines are ordered in the same way as the seeds
        seed_index = np.argmin(np.linalg.norm(line-np.expand_dims(seed, axis=0), axis=1))
        line[seed_index] = seed  # Some rounding happens on seed coordinates in the lines
        pairs.append([line[seed_index:], np.flip(line[:seed_index+1], axis=0)])
    return(pairs)

