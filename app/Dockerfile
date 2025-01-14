FROM python:3.12-slim-bullseye

ENV PYTHONUNBUFFERED 1

COPY . /app
WORKDIR /app
COPY ./scripts /scripts

EXPOSE 8000

RUN python3 -m venv /py

RUN apt-get update && apt-get install -y \
    build-essential \
    default-libmysqlclient-dev \
    && rm -rf /var/lib/apt/lists/*

ENV MYSQLCLIENT_CFLAGS="-I/usr/include/mysql"
ENV MYSQLCLIENT_LDFLAGS="-L/usr/lib/x86_64-linux-gnu -lmysqlclient"

# RUN build-base linux-headers && \ 
RUN /py/bin/pip install pip --upgrade && \
    /py/bin/pip install -r requirements.txt && \
    adduser --disabled-password --no-create-home app && \
    mkdir -p /vol/web/static && \   
    mkdir -p /vol/web/media && \
    chown -R app:app /vol && \
    chmod -R 755 /vol && \
    chmod -R +x /scripts

COPY ./static /vol/web/static
COPY ./media /vol/web/media


ENV PATH="/scripts:/py/bin:$PATH"

CMD ["gunicorn","proteoark.wsgi:application ","--bind", "0.0.0.0:8000"]
# CMD [ "run.sh" ]