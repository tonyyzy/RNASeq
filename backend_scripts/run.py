from configparser import ConfigParser
import os
from sqlalchemy import create_engine
from sqlalchemy.ext.automap import automap_base
from sqlalchemy.orm import Session
from uuid import UUID
import subprocess

if __name__ == "__main__":
    # change current working directory to project root
    root = os.path.realpath(os.path.join(os.path.dirname(__file__), '../..' ))
    os.chdir(root)

    # read config file
    config = ConfigParser()
    config.read("config.ini")
    database = config.get("main", "database")

    Base = automap_base()
    engine = create_engine(database)
    Base.prepare(engine, reflect=True)
    RSession = Base.classes.analysis_session
    Queue = Base.classes.analysis_queue
    session = Session(engine)

    query = session.query(RSession, Queue).filter(Queue.status == 1)\
                    .filter(RSession.id == Queue.session_id).all()
    if len(query) > 0:
        s, q = query[0]
        queue_id = q.id
        identifier = UUID(s.identifier)
        if q.jobtype == "index":
            logpath = f"{root}/Data/{identifier}/indexing.log"
            outdir = f"{root}/Data/{identifier}/genome"
        elif q.jobtype == "workflow":
            logpath = f"{root}/Data/{identifier}/workflow.log"
            outdir = f"{root}/Data/{identifier}/output"
        logfile = open(logpath, "a+")
        q.status = 0
        q.result = "submitted"
        session.commit()
        cwl = q.cwl
        yml = q.yml
        session.close()
        proc = subprocess.run(["cwl-runner",
                                f"--outdir={outdir}",
                                "--timestamp",
                                "--tmpdir-prefix=/tmp/",
                                "--tmp-outdir-prefix=/tmp/",
                                cwl,
                                yml],
                                stdout=logfile, stderr=logfile)

        Base = automap_base()
        engine = create_engine(database)
        Base.prepare(engine, reflect=True)
        Queue = Base.classes.analysis_queue
        session = Session(engine)
        q = session.query(Queue).filter(Queue.id == queue_id).first()
        if proc.returncode == 0:
            q.status = 0
            q.result = "success"
        else:
            q.status = 0
            q.result = "failed"
        
        session.commit()
        session.close()
        logfile.close()
