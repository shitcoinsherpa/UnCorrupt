"""Pytest fixtures for the corruption-detector tests."""
from __future__ import annotations

import os
import socket
import subprocess
import sys
import time
from collections.abc import Iterator

import pytest


def _wait_for_port(host: str, port: int, timeout: float = 60.0) -> None:
    deadline = time.time() + timeout
    while time.time() < deadline:
        try:
            with socket.create_connection((host, port), timeout=1):
                return
        except OSError:
            time.sleep(0.5)
    raise TimeoutError(f"server did not bind {host}:{port} within {timeout}s")


@pytest.fixture(scope="session")
def gradio_url() -> Iterator[str]:
    port = 7861
    env = {**os.environ, "GRADIO_SERVER_PORT": str(port)}
    proc = subprocess.Popen(
        [sys.executable, "-m", "uncorrupt.app"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        env=env,
    )
    try:
        _wait_for_port("127.0.0.1", port)
        yield f"http://127.0.0.1:{port}/"
    finally:
        proc.terminate()
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()
