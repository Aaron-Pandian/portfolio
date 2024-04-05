FROM python:3.10.12

COPY requirements.txt /app/requirements.txt
RUN pip install -r app/requirements.txt

COPY gene_api.py /app/gene_api.py
CMD ["python3", "app/gene_api.py"]