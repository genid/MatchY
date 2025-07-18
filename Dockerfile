# Use the official Python image
FROM python:3

# Set the working directory
WORKDIR /

# Copy only requirements first for caching
COPY requirements.txt .

# Install dependencies
RUN pip install --progress-bar off --no-cache-dir -r requirements.txt

# Copy the rest of the app
COPY . .

# Expose the Streamlit port
EXPOSE 8501

# Run the Streamlit app
CMD ["streamlit", "run", "1_🧬_Home.py", "--server.port=8501", "--server.address=0.0.0.0", "--server.fileWatcherType=none"]
