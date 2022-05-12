__version__ = "0.3.0"

from rich.console import Console

console = Console()

# from rich import inspect
# import snakemake

# console = Console()

# console.log(snakemake)
# inspect(snakemake)

# import time

# from rich.progress import Progress

# with Progress() as progress:

#     task1 = progress.add_task("[red]Downloading...", total=1000)
#     task2 = progress.add_task("[green]Processing...", total=1000)
#     task3 = progress.add_task("[cyan]Cooking...", total=1000)

#     while not progress.finished:
#         progress.update(task1, advance=0.5)
#         progress.update(task2, advance=0.3)
#         progress.update(task3, advance=0.9)
#         time.sleep(0.02)

# from time import sleep

# from rich.table import Column
# from rich.progress import Progress, BarColumn, TextColumn

# text_column = TextColumn("{task.description}", table_column=Column(ratio=1))
# bar_column = BarColumn(bar_width=None, table_column=Column(ratio=2))
# progress = Progress(text_column, bar_column, expand=True)

# with progress:
#     for n in progress.track(range(100)):
#         # progress.print(n)
#         sleep(0.1)
# from rich import box
# from rich.table import Table
# from rich.console import Console
# from snakemake.workflow import expand
# from pubsub import pub

# print(
#     expand(
#         "output/{sample}/fastqc_report/{sample}_{file_id}_fastqc.html",
#         sample=["A", "B"],
#         file_id=[1, 2],
#     )
# )

# console = Console()
# table = Table(title="Test table", box=box.ROUNDED, show_lines=True, expand=True)
# table.add_column("Job", justify="center", vertical="middle", no_wrap=True)
# table.add_column("Status", justify="center", no_wrap=True)
# sub_table = Table(
#     box=box.ROUNDED, show_lines=True, expand=True, show_header=False, show_edge=False
# )
# sub_table.add_column("Status", justify="center", no_wrap=True)
# sub_table.add_column("Time", justify="center", no_wrap=True)
# sub_table.add_row("Sample 1 started")
# sub_table.add_row("Sample 2 started")
# sub_table.add_row("Sample 3 finished", "5 min", style="black on green")
# # sub_table.add_row("Sample 3 finished [green]:heavy_check_mark:")
# table.add_row("Task1", sub_table)
# console.print(table)

# z = console.status("Monkeying around...", spinner="monkey")
# z.start()
# sleep(10)
# def listener1(arg1):
#     print("Function listener1 received:")
#     print("arg1 =", arg1)


# pub.subscribe(listener1, "rootTopic")

# print("Publish something via pubsub")
# anObj = dict(a=456, b="abc")
# pub.sendMessage("rootTopic", arg1=anObj)


"""Lite simulation of the top linux command."""
# import datetime
# import random
# import sys
# import time
# from dataclasses import dataclass

# from rich import box
# from rich.console import Console
# from rich.live import Live
# from rich.table import Table

# if sys.version_info >= (3, 8):
#     from typing import Literal
# else:
#     from typing_extensions import Literal


# @dataclass
# class Process:
#     pid: int
#     command: str
#     cpu_percent: float
#     memory: int
#     start_time: datetime.datetime
#     thread_count: int
#     state: Literal["running", "sleeping"]

#     @property
#     def memory_str(self) -> str:
#         if self.memory > 1e6:
#             return f"{int(self.memory/1e6)}M"
#         if self.memory > 1e3:
#             return f"{int(self.memory/1e3)}K"
#         return str(self.memory)

#     @property
#     def time_str(self) -> str:
#         return str(datetime.datetime.now() - self.start_time)


# def generate_process(pid: int) -> Process:
#     return Process(
#         pid=pid,
#         command=f"Process {pid}",
#         cpu_percent=random.random() * 20,
#         memory=random.randint(10, 200) ** 3,
#         start_time=datetime.datetime.now()
#         - datetime.timedelta(seconds=random.randint(0, 500) ** 2),
#         thread_count=random.randint(1, 32),
#         state="running" if random.randint(0, 10) < 8 else "sleeping",
#     )


# def create_process_table(height: int) -> Table:

#     processes = sorted(
#         [generate_process(pid) for pid in range(height)],
#         key=lambda p: p.cpu_percent,
#         reverse=True,
#     )
#     table = Table(
#         "PID", "Command", "CPU %", "Memory", "Time", "Thread #", "State", box=box.SIMPLE
#     )

#     for process in processes:
#         table.add_row(
#             str(process.pid),
#             process.command,
#             f"{process.cpu_percent:.1f}",
#             process.memory_str,
#             process.time_str,
#             str(process.thread_count),
#             process.state,
#         )

#     return table


# with Live(console=console, auto_refresh=False) as live:
#     while True:
#         live.update(create_process_table(console.size.height - 4), refresh=True)
#         time.sleep(1)

# import time

# from rich.progress import Progress

# with Progress() as progress:

#     task1 = progress.add_task("Downloading...", total=1000)
#     task2 = progress.add_task("Processing...", total=None)
#     task3 = progress.add_task("Cooking...", total=1000)

#     while not progress.finished:
#         progress.update(task1, advance=0.5)
#         progress.update(task2, advance=0.3)
#         progress.update(task3, advance=0.9)
#         time.sleep(0.02)


# from rich.table import Column

# text_column = TextColumn("{task.description}", table_column=Column(ratio=1))
# bar_column = BarColumn(bar_width=None, table_column=Column(ratio=2))
# progress = Progress(text_column, bar_column, expand=True)

# with progress:
#     for n in progress.track(range(100)):
#         progress.print(n)
#         sleep(0.1)

# with Progress() as progress:
#     task = progress.add_task("twiddling thumbs", total=10)
#     for job in range(10):
#         progress.console.print(f"Working on job #{job}")
#         sleep(2)
#         progress.advance(task)


# print(Panel("Hello, [red]World!", safe_box=True))


# import os
# from rich import print
# from time import sleep
# from rich.panel import Panel
# from rich.table import Table
# from rich.progress import Progress, BarColumn, TextColumn, TaskProgressColumn
# from rich.live import Live
# from rich.columns import Columns
# from rich.spinner import Spinner
# from rich.text import Text

# progress = Progress(
#     TextColumn("[progress.description]{task.description}"),
#     BarColumn(),
#     TaskProgressColumn(),
# )
# grid = Table.grid(expand=True)
# grid.add_column()
# for my_text in [
#     "Running FastQC [0/2 completed]",
#     "Running FastP [0/2 completed]",
#     "Running Trimmomatic [0/2 completed]",
# ]:
#     grid.add_row(Spinner("dots", style="blue", text=my_text))
# pipeline = progress.add_task("Running pipeline", total=100, completed=33)
# grid.add_row()
# grid.add_row(progress)

# live_display = Live(grid, refresh_per_second=20)
# live_display.start()

# sleep(2)
# grid = Table.grid(expand=True)
# grid.add_column()
# for my_text in [
#     "Running FastQC [1/2 completed]",
#     "Running FastP [1/2 completed]",
#     "Running Trimmomatic [1/2 completed]",
# ]:
#     grid.add_row(Spinner("dots", style="blue", text=my_text))
# progress.update(pipeline, total=100, completed=66)
# grid.add_row()
# grid.add_row(progress)
# live_display.update(grid)


# sleep(2)
# grid = Table.grid(expand=True)
# grid.add_column()
# for my_text in [
#     "✔ Completed FastQC [2/2 completed]",
#     "✔ Completed FastP [2/2 completed]",
#     "✔ Completed Trimmomatic [2/2 completed]",
# ]:
#     grid.add_row(Text(my_text, style="green"))
# progress.update(pipeline, total=100, completed=100)
# grid.add_row()
# grid.add_row(progress)
# live_display.update(grid)


# live_display.stop()

