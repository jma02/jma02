#!/usr/bin/env python3
"""
Script to generate programming jokes using the pyjokes library.
Includes error handling with a fallback joke.
"""

import pyjokes

def get_programming_joke():
    """
    Generate a programming joke using pyjokes library.
    Returns a fallback joke if pyjokes fails.
    """
    try:
        joke = pyjokes.get_joke()
        return joke
    except Exception as e:
        print(f"Error generating joke with pyjokes: {e}")
        # Fallback joke if pyjokes fails
        return "Why do programmers prefer dark mode? Because light attracts bugs!"

if __name__ == "__main__":
    joke = get_programming_joke()
    print(joke)