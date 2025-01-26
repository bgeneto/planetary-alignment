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
import concurrent.futures

BODY_MAP = {
    "mercury": "mercury",
    "venus": "venus",
    "mars": "mars",
    "jupiter": "jupiter",
    "saturn": "saturn",
    "uranus": "uranus",
    "neptune": "neptune",
    "moon": "moon",
}

TRANSLATIONS = {
    "en": {
        "title": "Planetary Alignment",
        "settings": "Settings",
        "select": "Select",
        "select_bodies": "Which objects do you want to check for alignment?",
        "alignment_method": "Alignment Detection Method",
        "pairwise_separation": "Maximum pairwise separation",
        "longitude_range": "Longitude range only",
        "alignment_threshold": "Alignment Threshold (degrees)",
        "max_pairwise_separation": "Max separation (degrees)",
        "location": "Observer Location",
        "latitude": "Latitude (degrees)",
        "longitude": "Longitude (degrees)",
        "date_range": "Date Range to Search",
        "start_date": "Start date",
        "end_year": "End year",
        "step_size_hours": "Time step (hours)",
        "check_visibility": "Check visibility at location",
        "min_altitude": "Minimum altitude (degrees)",
        "search_button": "Search for Alignment",
        "no_alignment_found": "No alignment found within the specified date range.",
        "alignment_found": "Alignment found: {}",
        "local_sky_positions": "Local Sky Positions (Alt-Az)",
        "bodies_warning": "Please select at least one celestial body to check for alignment.",
    },
    "pt": {
        "title": "Alinhamento Planetário",
        "settings": "Configurações",
        "select": "Selecione",
        "select_bodies": "Quais objetos você deseja verificar para alinhamento?",
        "alignment_method": "Método de detecção",
        "pairwise_separation": "Separação máxima entre pares",
        "longitude_range": "Apenas variação de longitude",
        "alignment_threshold": "Limite de alinhamento (graus)",
        "max_pairwise_separation": "Separação máxima (graus)",
        "location": "Localização do Observador",
        "latitude": "Latitude (graus)",
        "longitude": "Longitude (graus)",
        "date_range": "Faixa de Datas para Pesquisar",
        "start_date": "Data de início",
        "end_year": "Ano de término",
        "step_size_hours": "Intervalo temporal (horas)",
        "check_visibility": "Verificar visibilidade no local",
        "min_altitude": "Altura mínima (graus)",
        "search_button": "Pesquisar por Alinhamento",
        "no_alignment_found": "Nenhum alinhamento encontrado dentro da faixa de datas especificada.",
        "alignment_found": "Alinhamento encontrado: {}",
        "local_sky_positions": "Posições no Céu Local (Alt-Az)",
        "bodies_warning": "Nenhum corpo selecionado. Por favor, selecione pelo menos um.",
    },
}


def get_ecliptic_positions(bodies, date):
    time = Time(date)
    solar_system_ephemeris.set("de432s")
    ecliptic_positions = {}
    for body in bodies:
        coord = get_body(body, time)
        eclip_coord = coord.transform_to(GeocentricTrueEcliptic(equinox=time))
        ecliptic_positions[body] = {
            "lon": float(eclip_coord.lon.deg),
            "lat": float(eclip_coord.lat.deg),
        }
    return ecliptic_positions


def ecliptic_separation(lon1, lat1, lon2, lat2):
    lon1_rad, lat1_rad = np.radians([lon1, lat1])
    lon2_rad, lat2_rad = np.radians([lon2, lat2])
    cos_sep = np.sin(lat1_rad) * np.sin(lat2_rad) + np.cos(lat1_rad) * np.cos(
        lat2_rad
    ) * np.cos(lon1_rad - lon2_rad)
    cos_sep = np.clip(cos_sep, -1.0, 1.0)
    return np.degrees(np.arccos(cos_sep))


def calculate_ecliptic_separations(ecliptic_positions):
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


def check_ecliptic_alignment(positions, max_sep, method):
    if method == "pairwise":
        separations = calculate_ecliptic_separations(positions)
        return max(s[2] for s in separations) <= max_sep
    elif method == "longitude":
        longitudes = [pos["lon"] % 360 for pos in positions.values()]
        lon_range = max(longitudes) - min(longitudes)
        return min(lon_range, 360 - lon_range) <= max_sep


def get_horizontal_positions(bodies, date, location):
    time = Time(date)
    solar_system_ephemeris.set("de432s")
    altaz_positions = {}
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
    fig = go.Figure()
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
                r=[90 - pos["alt"]],
                theta=[pos["az"]],
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


def check_alignment_for_date(
    date, bodies, max_sep, method, check_vis, location, min_alt
):
    eclip_positions = get_ecliptic_positions(bodies, date)
    is_aligned = check_ecliptic_alignment(eclip_positions, max_sep, method)

    if check_vis and is_aligned:
        altaz = get_horizontal_positions(bodies, date, location)
        visible = all(pos["alt"] >= min_alt for pos in altaz.values())
        return (date, is_aligned and visible, eclip_positions)

    return (date, is_aligned, eclip_positions)


def main():
    lang_code = "en"
    st.set_page_config(
        page_title=TRANSLATIONS[lang_code]["title"], page_icon="💫", layout="wide"
    )

    lang = st.sidebar.selectbox("Language/Idioma", ["English", "Português"])
    lang_code = "en" if lang == "English" else "pt"

    st.title("💫" + TRANSLATIONS[lang_code]["title"])

    # --- Sidebar Controls ---
    st.sidebar.header(f":material/settings: {TRANSLATIONS[lang_code]['settings']}")

    st.sidebar.subheader(f"⭐️ {TRANSLATIONS[lang_code]['select_bodies']}")
    available_bodies = list(BODY_MAP.keys())
    selected_bodies = st.sidebar.multiselect(
        TRANSLATIONS[lang_code]["select"],
        available_bodies,
        default=["mercury", "venus", "mars", "jupiter", "saturn", "moon"],
    )
    if not selected_bodies:
        st.warning(TRANSLATIONS[lang_code]["bodies_warning"])
        return

    st.sidebar.subheader(f":telescope: {TRANSLATIONS[lang_code]['alignment_method']}")
    alignment_method = st.sidebar.radio(
        "",
        [
            TRANSLATIONS[lang_code]["pairwise_separation"],
            TRANSLATIONS[lang_code]["longitude_range"],
        ],
        index=0,
    )
    method = (
        "pairwise"
        if alignment_method == TRANSLATIONS[lang_code]["pairwise_separation"]
        else "longitude"
    )

    max_sep = st.sidebar.slider(
        TRANSLATIONS[lang_code]["alignment_threshold"],
        min_value=5.0,
        max_value=45.0,
        value=20.0,
        step=0.5,
    )

    st.sidebar.subheader(f":earth_americas: {TRANSLATIONS[lang_code]['location']}")
    lat = st.sidebar.number_input(
        TRANSLATIONS[lang_code]["latitude"], value=-15.7942, format="%0.4f"
    )
    lon = st.sidebar.number_input(
        TRANSLATIONS[lang_code]["longitude"], value=-47.8823, format="%0.4f"
    )
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)

    st.sidebar.subheader(f":calendar: {TRANSLATIONS[lang_code]['date_range']}")
    start_date = st.sidebar.date_input(
        TRANSLATIONS[lang_code]["start_date"], value=datetime.now().date()
    )
    end_year = st.sidebar.number_input(
        TRANSLATIONS[lang_code]["end_year"],
        min_value=start_date.year,
        max_value=start_date.year + 5,
        value=start_date.year + 1,
    )

    step_hours = st.sidebar.slider(
        TRANSLATIONS[lang_code]["step_size_hours"],
        min_value=1,
        max_value=24,
        value=6,
        help="Smaller steps increase accuracy but take longer to compute",
    )

    check_vis = st.sidebar.checkbox(
        TRANSLATIONS[lang_code]["check_visibility"], value=True
    )
    min_alt = (
        st.sidebar.slider(
            TRANSLATIONS[lang_code]["min_altitude"],
            min_value=-20,
            max_value=30,
            value=0,
            step=5,
        )
        if check_vis
        else 0
    )

    if st.button(TRANSLATIONS[lang_code]["search_button"]):
        with st.spinner("Searching..."):
            found_alignment = False
            all_dates = []
            current_date = datetime.combine(start_date, datetime.min.time())

            try:
                end_date = datetime(end_year, start_date.month, start_date.day)
            except ValueError:
                end_date = datetime(end_year, start_date.month, 28)

            while current_date <= end_date:
                all_dates.append(current_date)
                current_date += timedelta(hours=step_hours)

            api_bodies = [BODY_MAP[b] for b in selected_bodies]

            with concurrent.futures.ThreadPoolExecutor() as executor:
                futures = [
                    executor.submit(
                        check_alignment_for_date,
                        dt,
                        api_bodies,
                        max_sep,
                        method,
                        check_vis,
                        location,
                        min_alt,
                    )
                    for dt in all_dates
                ]

                for future in concurrent.futures.as_completed(futures):
                    date_checked, is_aligned, ecl_positions = future.result()
                    if is_aligned and not found_alignment:
                        found_alignment = True
                        st.success(
                            TRANSLATIONS[lang_code]["alignment_found"].format(
                                date_checked.strftime("%Y-%m-%d %H:%M")
                            )
                        )
                        altaz_positions = get_horizontal_positions(
                            api_bodies, date_checked, location
                        )
                        fig = create_sky_plot(altaz_positions)
                        st.plotly_chart(fig)
                        break

            if not found_alignment:
                st.error(TRANSLATIONS[lang_code]["no_alignment_found"])


if __name__ == "__main__":
    main()
