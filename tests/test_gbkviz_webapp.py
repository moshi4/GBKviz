import multiprocessing
import platform
import sys
import time

from gbkviz import gbkviz_webapp
from streamlit.cli import main as stcli_main


def test_gbkviz_webapp_launch():
    """test gbkviz webapp launch propery"""
    # TODO: Test only Linux now. This test method is not applicable to MacOS.
    #       Test MacOS as well in future.
    if platform.system() == "Linux":
        sys.argv = ["streamlit", "run", gbkviz_webapp.__file__]
        p = multiprocessing.Process(target=stcli_main, name="st")
        p.start()
        time.sleep(10)
        p.terminate()
        p.join()
        assert p.exitcode == 0
