.PHONY: setup pipeline dashboard

setup:
	python3 -m venv .venv
	.venv/bin/python -m pip install --upgrade pip
	.venv/bin/python -m pip install -r requirements.txt

pipeline:
	.venv/bin/python load_data.py
	PYTHONPATH=src MPLCONFIGDIR=.matplotlib XDG_CACHE_HOME=.cache .venv/bin/python -m teiko_analysis.pipeline

dashboard:
	PYTHONPATH=src .venv/bin/streamlit run dashboard/app.py
