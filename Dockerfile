FROM python:3.9-slim

ENV FLASK_APP=/usr/src/website/app.py
ENV PYTHONPATH=/usr/src/

COPY website/requirements.txt ./
RUN pip install --no-cache-dir --upgrade pip && \
      pip install --prefer-binary --no-cache-dir -r requirements.txt

COPY website/ /usr/src/website/
COPY database/ /usr/src/database/

WORKDIR /usr/src/website
CMD ["python3", "-m", "flask", "run", "--host=0.0.0.0", "--cert=adhoc", "--port", "8000"]
