import os
import signal
import logging
import threading
import flask.cli
from flask import Flask, request
from multiprocessing import Process
from .database import UpdateWorkflow, CreateDatabase, InitLiveDisplay, UpdateLiveDisplay

Server = None
Thread = threading.Lock()


def CreateServer(console, workflow, total):
    flask.cli.show_server_banner = lambda *args: None
    logging.getLogger("werkzeug").disabled = True

    app = Flask(__name__)
    # Getting the progress, pipeline and live display instance
    progress, pipeline, live_display = InitLiveDisplay(console, workflow, total)

    @app.route("/create", methods=["POST"])
    def Create():
        CreateDatabase(request.form.get("db_loc"))
        console.print()
        console.rule("[ Pipeline started running ]")
        # Live display started
        live_display.start()
        return "OK"

    @app.route("/print", methods=["POST"])
    def TerminalPrint():
        # Thread locking to prevent race condition
        Thread.acquire()
        # Updating the database
        update, total_completed = UpdateWorkflow(
            request.form.get("db_loc"),
            request.form.get("task"),
            request.form.get("sample"),
            request.form.get("status"),
            request.form.get("time"),
            workflow,
        )
        # Updating the live display
        UpdateLiveDisplay(
            update, total, progress, pipeline, live_display, total_completed
        )
        Thread.release()
        return "OK"

    @app.route("/exit", methods=["GET"])
    def ShutdownServer():
        # Stopping the live display
        live_display.stop()
        console.rule("[ Pipeline completed ]")
        # Killing the server itself
        os.kill(Server.pid, signal.SIGINT)
        return "OK"

    return app


def StartServer(host, port, debug, console, workflow, total):
    APP = CreateServer(console, workflow, total)
    APP.run(host=host, port=port, debug=debug)
    return


def StartLogger(console, workflow, total):
    global Server
    # Creating a separate process to handle the server and to get the Process ID for killing it
    Server = Process(
        target=StartServer, args=["127.0.0.1", 5000, False, console, workflow, total]
    )
    Server.start()
