import numpy as np
from astropy.time import Time
from astropy.coordinates import (
    get_body,
    EarthLocation,
    solar_system_ephemeris,
    SkyCoord,
    GeocentricTrueEcliptic,
)
import astropy.units as u
import streamlit as st
import plotly.graph_objects as go
from datetime import datetime, timedelta

import concurrent.futures  # Import for parallelization

# Mapping from user-friendly body names to the actual strings recognized by Astropy's get_body() function.
BODY_MAP = {
    "Mercury": "mercury",
    "Venus": "venus",
    "Mars": "mars",
    "Jupiter": "jupiter",
    "Saturn": "saturn",
    "Uranus": "uranus",
    "Neptune": "neptune",
    "Moon": "moon",
}

def get_ecliptic_positions(bodies, date):
    """
    For each entry in 'bodies', compute its geocentric ecliptic longitude and latitude
    at the specified date/time (a datetime.datetime object).

    Returns a dict: { body_name: {"lon": <float>, "lat": <float>} }
    """
    time = Time(date)
    # we can use a JPL ephemeris that includes outer planets
    solar_system_ephemeris.set('built-in')  # Em vez de 'de432s'

    ecliptic_positions = {}
    for body in bodies:
        # Get the sky coordinates for the body
        coord = get_body(body, time)
        # Transform to geocentric true ecliptic coordinates
        eclip_coord = coord.transform_to(GeocentricTrueEcliptic(equinox=time))
        ecliptic_positions[body] = {
            "lon": float(eclip_coord.lon.deg),
            "lat": float(eclip_coord.lat.deg),
        }

    return ecliptic_positions

def ecliptic_separation(lon1, lat1, lon2, lat2):
    """
    Compute spherical angular separation in degrees for two objects
    given their ecliptic coordinates (lon, lat in degrees).
    """
    # Convert to radians
    lon1_rad, lat1_rad = np.radians([lon1, lat1])
    lon2_rad, lat2_rad = np.radians([lon2, lat2])

    # Great-circle separation formula
    cos_sep = (
        np.sin(lat1_rad) * np.sin(lat2_rad)
        + np.cos(lat1_rad) * np.cos(lat2_rad) * np.cos(lon1_rad - lon2_rad)
    )
    cos_sep = np.clip(cos_sep, -1.0, 1.0)  # numerical safety
    return np.degrees(np.arccos(cos_sep))

def calculate_ecliptic_separations(ecliptic_positions):
    """
    Returns a list of pairwise (body1, body2, separation_degs)
    for all bodies in 'ecliptic_positions'.
    """
    bodies = list(ecliptic_positions.keys())
    pairs = []
    for i in range(len(bodies)):
        for j in range(i + 1, len(bodies)):
            b1, b2 = bodies[i], bodies[j]
            lon1, lat1 = ecliptic_positions[b1]["lon"], ecliptic_positions[b1]["lat"]
            lon2, lat2 = ecliptic_positions[b2]["lon"], ecliptic_positions[b2]["lat"]
            sep = ecliptic_separation(lon1, lat1, lon2, lat2)
            pairs.append((b1, b2, sep))
    return pairs

def check_ecliptic_alignment(positions, max_ecliptic_sep):
    """
    Checks if the maximum pairwise separation in ecliptic coordinates
    is <= max_ecliptic_sep.
    """
    separations = calculate_ecliptic_separations(positions)
    # The largest separation among all pairs:
    largest_sep = max(s[2] for s in separations)
    return largest_sep <= max_ecliptic_sep

def get_horizontal_positions(bodies, date, location):
    """
    Helper to create alt/az positions for plotting in a local sky plot.
    """
    time = Time(date)
    solar_system_ephemeris.set("de432s")

    altaz_positions = {}
    frame = EarthLocation(lat=location.lat, lon=location.lon).get_itrs(obstime=time)
    from astropy.coordinates import AltAz

    altaz_frame = AltAz(obstime=time, location=location)

    for body in bodies:
        coord = get_body(body, time)
        altaz = coord.transform_to(altaz_frame)
        altaz_positions[body] = {
            "alt": float(altaz.alt.deg),
            "az": float(altaz.az.deg),
        }
    return altaz_positions

def create_sky_plot(positions):
    """
    Creates a polar plot (Plotly) for alt-az positions of selected bodies.
    positions = { body_name: {"alt": float, "az": float} }
    """
    fig = go.Figure()
    # A simple color map if desired:
    color_map = {
        "mercury": "gray",
        "venus": "yellow",
        "mars": "red",
        "jupiter": "orange",
        "saturn": "gold",
        "uranus": "lightblue",
        "neptune": "blue",
        "moon": "white",
    }

    for body, pos in positions.items():
        label = body.capitalize()
        fig.add_trace(
            go.Scatterpolar(
                r=[90 - pos["alt"]],  # altitude -> distance from the zenith
                theta=[pos["az"]],  # azimuth
                mode="markers+text",
                text=[label],
                textposition="top center",
                name=label,
                marker=dict(size=10, color=color_map.get(body, "white")),
            )
        )

    fig.update_layout(
        polar=dict(
            radialaxis=dict(range=[0, 90], ticksuffix="°"),
            angularaxis=dict(direction="clockwise"),
        ),
        showlegend=True,
        title="Local Sky Positions (Alt-Az)",
    )

    return fig

# New helper function for parallel checking of alignment on a given date.
def check_alignment_for_date(date, bodies, max_sep):
    """
    Returns (date, is_aligned, eclip_positions).
    """
    eclip_positions = get_ecliptic_positions(bodies, date)
    is_aligned = check_ecliptic_alignment(eclip_positions, max_sep)
    return (date, is_aligned, eclip_positions)

def main():
    st.set_page_config(page_title="Planetary Alignment", page_icon="💫", layout="wide")
    st.title("💫 Planetary (and Moon) Alignment")

    # --- Sidebar Controls ---
    st.sidebar.header(":material/settings: Configuration")

    # 1) Select which bodies to consider:
    st.sidebar.subheader("💫 Which objects do you want to check for alignment?")
    available_bodies = list(BODY_MAP.keys())
    selected_bodies = st.sidebar.multiselect(
        "Select bodies",
        available_bodies,
        default=[
            "Mercury",
            "Venus",
            "Mars",
            "Jupiter",
            "Saturn",
            "Moon",
        ],  # Default selection
    )
    if not selected_bodies:
        st.warning("No bodies selected. Please select at least one.")
        return

    # 2) Ecliptic alignment threshold:
    st.sidebar.subheader(":telescope: Alignment Threshold (Ecliptic)")
    max_ecliptic_sep = st.sidebar.slider(
        "Max pairwise separation (degrees)",
        min_value=0.0,
        max_value=30.0,
        value=10.0,
        step=0.5,
        help="Maximum angular separation in ecliptic coords to consider 'aligned'.",
    )

    # 3) Location
    st.sidebar.subheader(":earth_americas: Observer Location")
    lat = st.sidebar.number_input("Latitude (degrees)", value=-15.7942, format="%0.4f")
    lon = st.sidebar.number_input("Longitude (degrees)", value=-47.8823, format="%0.4f")
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)

    # 4) Date range
    st.sidebar.subheader(":calendar: Date Range to Search")
    start_date = st.sidebar.date_input("Start date", value=datetime.now().date())
    end_year = st.sidebar.number_input(
        "End year", min_value=start_date.year, max_value=2100, value=start_date.year + 1
    )

    # We fix the month/day to be the same as start_date for the end date
    try:
        end_datetime = datetime(end_year, start_date.month, start_date.day)
    except ValueError:
        end_datetime = datetime(end_year, start_date.month, min(start_date.day, 28))

    start_datetime = datetime.combine(start_date, datetime.min.time())
    step_days = st.sidebar.slider("Step size (days)", 1, 30, 1)

    # --- Search Button ---
    if st.button("Search for Alignment"):
        with st.spinner("Searching..."):
            found_alignment = False
            # Create a list of all the dates we will check
            all_dates = []
            current_date = start_datetime
            while current_date <= end_datetime:
                all_dates.append(current_date)
                current_date += timedelta(days=step_days)

            # Convert user-chosen bodies to the actual strings used by get_body
            api_bodies = [BODY_MAP[b] for b in selected_bodies]

            # Parallel check each date for alignment
            with concurrent.futures.ThreadPoolExecutor() as executor:
                results = list(
                    executor.map(
                        lambda dt: check_alignment_for_date(dt, api_bodies, max_ecliptic_sep),
                        all_dates
                    )
                )

            # Now find the earliest date that is aligned
            for (date_checked, is_aligned, ecl_positions) in results:
                if is_aligned:
                    found_alignment = True
                    st.success(f"Alignment found: {date_checked.strftime('%Y-%m-%d %H:%M')}")
                    # For visualization: get alt/az on that date
                    altaz_positions = get_horizontal_positions(api_bodies, date_checked, location)
                    fig = create_sky_plot(altaz_positions)
                    st.plotly_chart(fig)
                    break  # Stop as soon as we find the first alignment

            if not found_alignment:
                st.error("No alignment found within the specified date range.")

if __name__ == "__main__":
    main()
