FROM python:3.12-slim

WORKDIR /app
COPY pyproject.toml README.md ./
COPY src ./src
RUN pip install --no-cache-dir .

ENV MOLBIO_DATA_DIR=/var/lib/localmolbio
VOLUME ["/var/lib/localmolbio"]
EXPOSE 8000
CMD ["uvicorn", "localmolbio.api:app", "--host", "0.0.0.0", "--port", "8000"]
