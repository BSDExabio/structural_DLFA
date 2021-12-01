FROM python:3.8-slim

ENV FLASK_APP=/usr/src/app.py

WORKDIR /usr/src/

COPY requirements.txt ./
RUN pip install --no-cache-dir --upgrade pip && \
      pip install --prefer-binary --no-cache-dir -r requirements.txt

COPY . .

CMD ["flask", "run", "--host=0.0.0.0"]
