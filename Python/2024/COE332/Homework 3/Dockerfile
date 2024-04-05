FROM ubuntu:20.04

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y python3 && \
    apt-get install -y python3-pip

RUN pip3 install pytest==8.0.0

RUN pip3 install matplotlib && \
    pip3 install scipy && \
    pip3 install pytest

COPY ml_data_analysis.py /code/ml_data_analysis.py
COPY test_ml_data_analysis.py /code/test_ml_data_analysis.py
COPY gcd_algorithm.py /code/gcd_algorithm.py
COPY test_gcd_algorithm.py /code/test_gcd_algorithm.py

RUN chmod +rx /code/ml_data_analysis.py && \
    chmod +rx /code/test_ml_data_analysis.py && \
    chmod +rx /code/test_gcd_algorithm.py && \
    chmod +rx /code/gcd_algorithm.py

ENV PATH="/code:$PATH"
