#!/usr/bin/env python3
"""
Script to generate a random joke using the pyjokes library.
This script will be executed by the GitHub Actions workflow.
"""

import pyjokes


def generate_joke():
    """Generate a random joke using pyjokes library."""
    try:
        joke = pyjokes.get_joke()
        return joke
    except Exception as e:
        # Fallback joke in case of any issues
        return "Why do programmers prefer dark mode? Because light attracts bugs!"


if __name__ == "__main__":
    joke = generate_joke()
    print(joke)