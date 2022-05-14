import sqlite3
from rich import box
from time import sleep
from pathlib import Path
from sqlite3 import Error
from rich.live import Live
from rich.text import Text
from rich.table import Table
from datetime import datetime
from rich.spinner import Spinner
from rich.console import Console
from rich.progress import (
    Progress,
    BarColumn,
    TextColumn,
    TaskProgressColumn,
    TimeElapsedColumn,
)


def CreateConnection(db_loc):
    db_file = "db.sqlite3"
    connection = None
    try:
        connection = sqlite3.connect(Path(db_loc, db_file))
    except Error as e:
        print(e)
    finally:
        return connection


def CreateDatabase(db_loc):
    """
    Create a database connection to a SQLite database
    """
    connection = CreateConnection(db_loc)
    if connection:
        connection.execute(
            """
                CREATE TABLE "NITRO"
                    (
                        "id" integer NOT NULL PRIMARY KEY AUTOINCREMENT,
                        "Job" varchar(255) NOT NULL,
                        "Sample" varchar(255) NOT NULL,
                        "Status" varchar(255) NOT NULL,
                        "Time" datetime NOT NULL
                    );
            """
        )
        connection.close()


def UpdateDatabase(db_loc, job, sample, status, time):
    connection = CreateConnection(db_loc)
    if connection:
        connection.execute(
            f"""
                INSERT INTO "NITRO" ("Job", "Sample", "Status", "Time")
                    VALUES (?,?,?,?)
            """,
            (job, sample, status, time),
        )
        connection.commit()
        connection.close()


def ListDatabase(db_loc):
    connection = CreateConnection(db_loc)
    if connection:
        connection.row_factory = lambda cursor, row: {
            "Job": row[1],
            "Sample": row[2],
            "Status": row[3],
            "Time": datetime.strptime(row[4], "%Y-%m-%d %H:%M:%S.%f"),
        }
        crsr = connection.cursor()
        crsr.execute("SELECT * FROM 'NITRO'")
        result = crsr.fetchall()
        response = {row["Job"]: {} for row in result}
        for row in result:
            response[row["Job"]][row["Sample"]] = row["Status"]
        connection.close()
        return response


def UpdateWorkflow(db_loc, job, sample, status, time, workflow):
    UpdateDatabase(db_loc, job, sample, status, time)
    progress = ListDatabase(db_loc)
    total_completed = 0
    order = []
    for tool, value in progress.items():
        if not tool in order:
            order.append(tool)
        running = len([item for item, state in value.items() if (state == "Started")])
        completed = len(
            [item for item, state in value.items() if (state == "Finished")]
        )
        workflow[tool]["running"] = running
        workflow[tool]["completed"] = completed
        total_completed = total_completed + completed

    # Ordering the workflow
    new_workflow = {}
    # Order from the database
    for item in order:
        new_workflow[item] = workflow[item]
    # Put rest of the items in any order
    for tool, value in workflow.items():
        if not tool in list(new_workflow.keys()):
            new_workflow[tool] = value
    return new_workflow, total_completed


def InitLiveDisplay(console, workflow, total):
    """
    Setting the progress bar style with
        - Text column (message)
        - Bar column (progress bar)
        - Task progress column (percentage)
        - Time elapsed column (time elapsed)
    """
    progress = Progress(
        TextColumn("[progress.description]{task.description}"),
        BarColumn(bar_width=int(console.size.width * 0.9)),
        TaskProgressColumn(),
        TimeElapsedColumn(),
    )
    grid = Table.grid(expand=True)
    grid.add_column()
    grid.add_row()
    # Adding text for not started tools
    for Tool, Status in workflow.items():
        grid.add_row(Text(f"• {Tool}: Not started", style="dim italic"))
    # Creating the progress bar instance and saving it
    pipeline = progress.add_task("Running pipeline", total=total, completed=0)
    grid.add_row()
    grid.add_row(progress)
    grid.add_row()
    # Creating a live display and aving its instance
    live_display = Live(grid, console=console, refresh_per_second=60)

    return progress, pipeline, live_display


def UpdateLiveDisplay(
    workflow, total, progress, pipeline, live_display, total_completed
):
    grid = Table.grid(expand=True)
    grid.add_column()
    grid.add_row()
    for Tool, Status in workflow.items():
        # Checking whether the tool has completed the work
        if Status["completed"] == Status["total"]:
            text = f"✔ {Tool}: Completed"
            grid.add_row(Text(text, style="green"))
        # Checking whether tool has not started
        elif Status["completed"] == 0 and Status["running"] == 0:
            grid.add_row(Text(f"• {Tool}: Not started", style="dim italic"))
        # Tool has started and has it has either some runnning instance or completed instance
        else:
            text = (
                f"[yellow][italic]{Tool}[/italic][/yellow]: "
                + f"[dim]({Status['running']} Running)[/dim] "
                + f"[bold green]({Status['completed']}/{Status['total']} Completed)[/bold green]"
            )
            grid.add_row(Spinner("dots", style="blue", text=text))
    progress.update(pipeline, total=(total - 1), completed=total_completed)
    grid.add_row()
    grid.add_row(progress)
    grid.add_row()
    # Updating the live display
    live_display.update(grid)
