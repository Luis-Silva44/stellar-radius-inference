from astroquery.utils.tap.core import TapPlus


SWEETCAT_URL = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap"
SWEETCAT_TABLE = '"J/A+A/691/A53/catalog"'


def query_sweetcat(star_name):
    tap = TapPlus(url=SWEETCAT_URL)

    query = f"""
        SELECT
            "Name",
            "Teff",
            "e_Teff",
            "logg",
            "e_logg",
            "[Fe/H]",
            "e_[Fe/H]",
            "Radius-t",
            "e_Radius-t",
            "GaiaDR3",
            "Plx",
            "e_Plx",
            "Dist"
        FROM {SWEETCAT_TABLE}
        WHERE "Name" = '{star_name}'
    """

    job = tap.launch_job(query)
    return job.get_results()


