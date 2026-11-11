@echo off
REM Double-click this file to start YTAB. It needs Docker Desktop and nothing else: R, Shiny, Python and the
REM aligners all live inside the image, and your data never leaves this computer.
REM
REM macOS and Linux: use YTAB.command instead.

setlocal
set IMAGE=ghcr.io/jayelazuno/ytab-app:v0.1
cd /d "%~dp0"
set REPO=%CD%

echo YTAB - Yeast Transposon Analysis Browser
echo ----------------------------------------

if not exist "%REPO%\app\shiny\run_app.R" (
  echo.
  echo ERROR: This file must stay inside the YTAB folder; it is looking in %REPO%.
  goto :halt
)

where docker >nul 2>&1
if errorlevel 1 (
  echo.
  echo ERROR: Docker is not installed.
  echo Install Docker Desktop from https://docs.docker.com/get-started/get-docker/ , start it,
  echo then double-click this again.
  goto :halt
)

docker info >nul 2>&1
if errorlevel 1 (
  echo.
  echo ERROR: Docker is installed but not running.
  echo Start Docker Desktop, wait for it to say it is running, then double-click this again.
  goto :halt
)

docker image inspect %IMAGE% >nul 2>&1
if errorlevel 1 (
  echo.
  echo First run: downloading YTAB ^(about 1.4 GB^). This happens once and takes a few minutes...
  docker pull %IMAGE%
  if errorlevel 1 (
    echo.
    echo ERROR: The download failed. Check your internet connection and try again.
    goto :halt
  )
)

set PORT=3838
set NAME=ytab-%PORT%
docker rm -f %NAME% >nul 2>&1

echo.
echo Starting YTAB on http://localhost:%PORT%
REM Your user folder is handed to YTAB read-only at /host, so "Start from FASTQ" can browse to your own
REM sequencing files. YTAB only ever reads them; everything it writes goes to this YTAB folder.
docker run -d --name %NAME% -p %PORT%:3838 -v "%REPO%:/opt/ytab" -v "%USERPROFILE%:/host:ro" %IMAGE% >nul
if errorlevel 1 (
  echo.
  echo ERROR: YTAB could not start. Run this for the reason:  docker logs %NAME%
  goto :halt
)

echo Waiting for it to be ready...
timeout /t 25 /nobreak >nul
start "" "http://localhost:%PORT%"

echo.
echo YTAB is running at http://localhost:%PORT%
echo Your results are saved in this folder, under output\.
echo.
echo To analyse your own reads, use "Start from FASTQ" and click "Browse local folders":
echo your files are there under "My files", read-only.
echo.
echo Leave this window open while you work.
echo To stop YTAB, close this window or press Control-C.
echo.
docker logs -f %NAME%

docker rm -f %NAME% >nul 2>&1
goto :eof

:halt
echo.
pause
