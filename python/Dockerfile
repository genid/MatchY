FROM python:3.12-slim-bullseye

WORKDIR /

RUN apt-get update && apt-get install -y \
    build-essential \
    graphviz \
    graphviz-dev \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install --progress-bar off --no-cache-dir -r requirements.txt

COPY . .

EXPOSE 8501
CMD ["streamlit", "run", "1_🧬_Home.py", "--server.port=8501", "--server.address=0.0.0.0", "--server.fileWatcherType=none"]
