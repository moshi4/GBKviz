import multiprocessing
import sys
import time

from gbkviz import gbkviz_webapp
from streamlit.cli import main as stcli_main


def test_gbkviz_webapp_launch():
    """test gbkviz webapp launch propery"""
    sys.argv = ["streamlit", "run", gbkviz_webapp.__file__]
    p = multiprocessing.Process(target=stcli_main, name="st")
    p.start()
    time.sleep(5)
    p.terminate()
    p.join()
    assert p.exitcode == 0
