import pandas as pd
import gpxpy
import folium
import numpy as np
import geopy.distance
import pyproj
import csv

def get_data_from_gpx(file):
    gpx_file = open(file, "r")
    gpx = gpxpy.parse(gpx_file)
    coordinates = []
    for track in gpx.tracks:
        for segment in track.segments:
            for point in segment.points:
                coordinates.append([point.latitude, point.longitude, point.elevation, point.time, -1])

    coordinates[0].append(-9999999)
    coordinates[0].append(-9999999)
    for i in range(1, len(coordinates)):
        coordinates[i].append(geopy.distance.geodesic(coordinates[i][0:2], coordinates[i -1][0:2]).km)
        time_dif_delta = coordinates[i][3] - coordinates[i - 1][3]
        time_dif = time_dif_delta.total_seconds()
        coordinates[i].append((coordinates[i][5] / time_dif)*60*60) #km/h

    pts_np = np.array(coordinates)

    return pts_np

def mark_uphill(pts):
    uphill_starts = []
    uphill_ends = []
    is_lift = False
    downhill_count = 0
    for i in range(1, len(pts)):
        if pts[i][1] == -111.551275:
            breakpoint
        if pts[i][2] >= pts[i-1][2]:
            is_lift = True
            pts[i-1][4] = 1
            uphill_starts.append([pts[i - 1][0], pts[i - 1][1]])
        elif pts[i][2] + 0.5 < pts[i-1][2]:
            downhill_count += 1
            if downhill_count < 2 and is_lift:
                pts[i-1][4] = 1
            else:
                downhill_count = 0
                is_lift = False
                pts[i-1][4] = 0
                uphill_ends.append([pts[i - 1][0], pts[i - 1][1]])
        else:
            if is_lift:
                pts[i-1][4] = 1
            else:
                pts[i-1][4] = 0

        if i == len(pts) - 1:
            pts[i][4] = pts[i -1][4]
                

    return pts, uphill_starts, uphill_ends

def unmark_small_uphills(pts):
    i = 0
    uphill_start_index = None
    is_uphill = False
    while i < len(pts) - 1:
        if pts[i][4] == 1:
            if uphill_start_index is None:
                uphill_start_index = i
            is_uphill = True
        elif is_uphill and pts[i][4] == 0:
            uphill_dist = geopy.distance.geodesic(pts[i - 1][0:2], pts[uphill_start_index][0:2]).km
            # Marks uphill segment as not a lift if shorter than 0.12km
            if uphill_dist < 0.1:
                for j in range(uphill_start_index, i):
                    pts[j][4] = 0
            uphill_start_index = None
            is_uphill = False
        i += 1

    return pts


def mark_small_downhills(pts):
    i = 0
    downhill_start_index = None
    is_lift = True
    while i < len(pts) - 1:
        if pts[i][4] == 1:
            if is_lift is False:
                downhill_dist = geopy.distance.geodesic(pts[i - 1][0:2], pts[downhill_start_index][0:2]).km
                if downhill_dist < 0.2:
                    for j in range(downhill_start_index, i):
                        pts[j][4] = 1
            is_lift = True
        elif pts[i][4] == 0 and is_lift is True:
            downhill_start_index = i
            is_lift = False

        i += 1

    return pts


def mark_downhill_lifts(pts):
    geodesic = pyproj.Geod(ellps='WGS84')
    i = 0
    uphill_start_index = None
    uphill_end_index = None
    is_uphill = False
    bearing_list = []
    while i < len(pts) - 1:
        if pts[i][4] == 1:
            if uphill_start_index is None:                
                if uphill_end_index is not None:
                    # TODO: Add more points for checking bearings. (Snowbird ^Little Cloud vMineral counts as lift)
                    down_fwd_bearing, down_bwd_bearing, down_dist = geodesic.inv(pts[uphill_end_index][1], pts[uphill_end_index][0], pts[i][1], pts[i][0])
                    bearing_dif = down_fwd_bearing - up_fwd_bearing
                    if abs(bearing_dif) < 10:
                        bearing_list.append([up_fwd_bearing, up_bwd_bearing, down_fwd_bearing, down_bwd_bearing, pts[uphill_end_index][1], pts[uphill_end_index][0], pts[i][1], pts[i][0]])
                        for j in range(uphill_end_index + 1, i):
                            pts[j][4] = 1

                uphill_start_index = i
            uphill_end_index = None
            is_uphill = True
        elif is_uphill and pts[i][4] == 0:
            if uphill_end_index is None:
                if uphill_start_index is not None:
                    up_fwd_bearing, up_bwd_bearing, up_dist = geodesic.inv(pts[uphill_start_index][1], pts[uphill_start_index][0], pts[i-1][1], pts[i-1][0])
                uphill_end_index = i - 1

            uphill_start_index = None
            is_uphill = False

        i += 1

    with open("out_bearing.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(bearing_list)

    return pts, bearing_list


# fwd_azimuth,back_azimuth,distance = geodesic.inv(long1, lat1, long2, lat2)


def generate_map_multiple(bearings):
    ski_map = folium.Map()

    vector_coords = []
    for i in range(0, len(bearings)):
        vector_coords.append([[bearings[i][5], bearings[i][4]],[bearings[i][7],bearings[i][6]]])

    for vect in vector_coords:
        folium.PolyLine(vect, weight=2, opacity=1).add_to(ski_map)
        folium.Marker(vect[0], icon=folium.Icon(color="green")).add_to(ski_map)

    ski_map.save("test_0.html")


def generate_map(new_coords, old_coords, pts):
    ski_map = folium.Map()
    folium.PolyLine(old_coords, weight=3, opacity=1, color="red").add_to(ski_map)

    folium.PolyLine(new_coords, weight=2, opacity=1).add_to(ski_map)
    for pt in pts:
        folium.Marker([pt[0], pt[1]], icon=folium.Icon(color="green"), popup=pt).add_to(ski_map)
    ski_map.save("test_1.html")

def get_non_lift_coordinates(pts):
    pts_no_lift = []
    for i in range(0, len(pts)):
        if pts[i][4] == 0:
            pts_no_lift.append(pts[i][0:2])

    return pts_no_lift

def save_to_csv(pts):
    with open("out_data.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(pts)


def get_stop_start_pts(pts):
    new_pts = []
    is_lift = False
    for pt in pts:
        if is_lift and pt[4] == 0:
            is_lift = False
            new_pts.append(pt)
        elif not is_lift and pt[4] == 1:
            is_lift = True
            new_pts.append(pt)

    return new_pts


def main():
    gpx_file_name = "Day_10_The_Bird.gpx"
    gpx_data = get_data_from_gpx(gpx_file_name)
    old_coordinates = gpx_data[:,[0,1]]
    gpx_data, start_pts, end_pts = mark_uphill(gpx_data)
    gpx_data = unmark_small_uphills(gpx_data)
    gpx_data, bearings = mark_downhill_lifts(gpx_data)
    gpx_data = mark_small_downhills(gpx_data)
    save_to_csv(gpx_data)
    coordinates_no_lift = get_non_lift_coordinates(gpx_data)
    ss_data = get_stop_start_pts(gpx_data)
    generate_map(coordinates_no_lift, old_coordinates, ss_data)
    # new_coordinates = new_points[:,[0,1]]
    # generate_map_multiple(bearings)
    breakpoint



if __name__ == "__main__":
    main()