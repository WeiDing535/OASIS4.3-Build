"""OASIS4.3 focused local application."""
def main():
    from .gui import main as serve
    serve(focused=True)

if __name__=='__main__':main()
