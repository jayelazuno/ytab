#!/usr/bin/env python3
import argparse,socket,subprocess,sys,time,urllib.request
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
def free_port(host):
 with socket.socket() as s:s.bind((host,0));return s.getsockname()[1]
def main():
 p=argparse.ArgumentParser(description="Start the local Shiny app, verify HTTP, and stop it cleanly.");p.add_argument("--project-config",type=Path,required=True);p.add_argument("--host",default="127.0.0.1");p.add_argument("--port",type=int);p.add_argument("--timeout",type=float,default=45);a=p.parse_args();port=a.port or free_port(a.host);cfg=a.project_config.resolve();url=f"http://{a.host}:{port}";cmd=["Rscript",str(ROOT/"app/shiny/run_app.R"),"--host",a.host,"--port",str(port),"--project-config",str(cfg),"--no-browser"];start=time.monotonic();proc=subprocess.Popen(cmd,cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
 try:
  while time.monotonic()-start<a.timeout:
   if proc.poll() is not None:raise RuntimeError("Shiny exited during startup")
   try:
    with urllib.request.urlopen(url,timeout=1) as response:
     html=response.read().decode("utf-8",errors="replace")
     required=("YTAB — Yeast Transposon Analysis Browser",'id="app_shell"')
     if response.status<500 and all(text in html for text in required):print(f"PASS\nLanding-page shell found: {', '.join(required)}\nURL tested: {url}\nStartup elapsed seconds: {time.monotonic()-start:.3f}");return 0
   except Exception:time.sleep(.25)
  raise TimeoutError(f"No HTTP response within {a.timeout} seconds")
 except Exception as e:
  print(f"FAIL: {e}\nURL tested: {url}",file=sys.stderr);return 1
 finally:
  if proc.poll() is None:proc.terminate()
  try:out,err=proc.communicate(timeout=8)
  except subprocess.TimeoutExpired:proc.kill();out,err=proc.communicate()
  if proc.returncode not in (0,-15) and err:print("Captured stderr:\n"+err,file=sys.stderr)
if __name__=="__main__":raise SystemExit(main())
