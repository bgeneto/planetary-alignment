# Planetary Alignment

A Python script to search for alignments of planets and the Moon in the ecliptic plane.

## Description

This script uses the Astropy library to calculate the positions of planets and the Moon in the ecliptic plane and checks for alignments within a specified date range. The alignment threshold can be adjusted, and the script also generates a local sky plot (alt-az) for visualization.

## Features

* Search for alignments of planets and the Moon in the ecliptic plane
* Adjustable alignment threshold (max pairwise separation in ecliptic coordinates)
* Local sky plot (alt-az) generation for visualization
* User-friendly interface using Streamlit

## Requirements

* Python 3.x
* Astropy library
* Streamlit library
* Plotly library

## Usage

1. Clone the repository
2. Install the required libraries using `pip install -r requirements.txt`
3. Run the script using `python main.py`
4. Follow the prompts in the Streamlit interface to configure the search parameters

## Example Use Case

* Search for alignments of Mercury, Venus, and Mars within a 10-degree threshold in the ecliptic plane between 2022-01-01 and 2023-01-01.

## Contributing

Contributions are welcome! If you'd like to add new features or improve the existing code, please submit a pull request.
