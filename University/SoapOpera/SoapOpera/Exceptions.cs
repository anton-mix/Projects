using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SoapOpera {
    public class OutOfBorderPhotoException : Exception {
        public string exception = "";
        public OutOfBorderPhotoException() {}

        public OutOfBorderPhotoException(string exception) {
            this.exception = exception;
        }
    }

    public class OutOfBorderCommentsException : Exception {
        public string exception = "";
        public OutOfBorderCommentsException() {}

        public OutOfBorderCommentsException(string exception) {
            this.exception = exception;
        }
    }

    public class ErrorFindPhotoException : Exception {
        public string TitleAlbum, exception;
        public ErrorFindPhotoException(string TitleAlbum) {
            this.TitleAlbum = TitleAlbum;
        }
        public ErrorFindPhotoException(string TitleAlbum, string exception) {
            this.TitleAlbum = TitleAlbum;
            this.exception = exception;
        }
    }

    public class ErrorFindCommentsException : Exception {
        public int NumberPhoto;
        public string TitleAlbum, exception;
        public ErrorFindCommentsException(int NumberPhoto, string TitleAlbum) {
            this.NumberPhoto = NumberPhoto;
            this.TitleAlbum = TitleAlbum;
        }

        public ErrorFindCommentsException(int NumberPhoto, string TitleAlbum, string exception) {
            this.NumberPhoto = NumberPhoto;
            this.TitleAlbum = TitleAlbum;
            this.exception = exception;
        }
    }

    public class ErrorFindNameUserException : Exception {
        public string TitleAlbum, exception;
        public ErrorFindNameUserException(string TitleAlbum) {
            this.TitleAlbum = TitleAlbum;
        }
        public ErrorFindNameUserException(string TitleAlbum, string exception) {
            this.TitleAlbum = TitleAlbum;
            this.exception = exception;
        }
    }

    public class NotCommentException : Exception {
        public string exception = "";
        public NotCommentException() { }

        public NotCommentException(string exception) {
            this.exception = exception;
        }
    }

    public class BadFormatException : Exception {
        public string Error, exception;
        public BadFormatException(string Error) {
            this.Error = Error;
        }

        public BadFormatException(string Error, string exception) {
            this.Error = Error;
            this.exception = exception;
        }
    }

    public class WarrningFormatId : Exception {
        public string exception = "";
        public WarrningFormatId() { }

        public WarrningFormatId(string exception) {
            this.exception = exception;
        }
    }

    public class WarrningFormatCount : Exception {
        public string exception = "";
        public WarrningFormatCount() { }

        public WarrningFormatCount(string exception) {
            this.exception = exception;
        }
    }

    public class WarrningFormatName : Exception {
        public string exception = "";
        public WarrningFormatName() { }

        public WarrningFormatName(string exception) {
            this.exception = exception;
        }
    }

    public class UnknownColumnException : Exception {
        public string exception = "";
        public UnknownColumnException() { }

        public UnknownColumnException(string exception) {
            this.exception = exception;
        }
    }

    public class BadFormatSafeException : Exception {
        public string exception = "";
        public BadFormatSafeException() { }

        public BadFormatSafeException(string exception) {
            this.exception = exception;
        }
    }

}
