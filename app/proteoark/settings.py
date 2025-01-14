import os
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent

SECRET_KEY = "django-insecure-xji_t-5+)%y!fs#odgz83dkg@la4ee!u36qs$hf557$6qdblc7"

DEBUG = True

ALLOWED_HOSTS = [
    "ciods.in",
    "proteoark.ciods.in",
    "localhost",
    "127.0.0.1",
    "app",        # Docker service name
]
# Application definition

CSRF_TRUSTED_ORIGINS = ['http://localhost:3000', 
                         "https://ciods.in",  
                            "https://proteoark.ciods.in",
                         'http://localhost:8000','http://localhost:80']

CORS_EXPOSE_HEADERS=['Content-Type', 'X-CSRFToken']
SESSION_COOKIE_SECURE = True
SESSION_COOKIE_SAMESITE = 'Lax'
CSRF_COOKIE_SAMESITE = 'Lax'
CSRF_COOKIE_SECURE = True
CSRF_COOKIE_HTTPONLY = True

INSTALLED_APPS = [
    'whitenoise.runserver_nostatic',

    "django.contrib.admin",
    "django.contrib.auth",
    "django.contrib.contenttypes",
    "django.contrib.sessions",
    "django.contrib.messages",
    "django.contrib.staticfiles",
    
    'proteome',
]

MIDDLEWARE = [
    'django.middleware.security.SecurityMiddleware',
    'whitenoise.middleware.WhiteNoiseMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'corsheaders.middleware.CorsMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
]

ROOT_URLCONF = "proteoark.urls"

TEMPLATES = [
    {
        "BACKEND": "django.template.backends.django.DjangoTemplates",
        'DIRS': [BASE_DIR / 'templates'],
        "APP_DIRS": True,
        "OPTIONS": {
            "context_processors": [
                "django.template.context_processors.debug",
                "django.template.context_processors.request",
                "django.contrib.auth.context_processors.auth",
                "django.contrib.messages.context_processors.messages",
            ],
        },
    },
]

CORS_ORIGIN_WHITELIST = [
    'http://localhost:3000',  # Replace with your Next.js frontend URL
    "https://ciods.in",
    "https://proteoark.ciods.in",
    'http://localhost:8000',
]

CORS_ALLOW_CREDENTIALS = True
CORS_ORIGIN_ALLOW_ALL = True

WSGI_APPLICATION = "proteoark.wsgi.application"

DATABASES = {
    "default": {
        "ENGINE": "django.db.backends.sqlite3",
        "NAME": BASE_DIR / "db.sqlite3",
    }
}


# Password validation
# https://docs.djangoproject.com/en/5.0/ref/settings/#auth-password-validators

AUTH_PASSWORD_VALIDATORS = [
    {
        "NAME": "django.contrib.auth.password_validation.UserAttributeSimilarityValidator",
    },
    {
        "NAME": "django.contrib.auth.password_validation.MinimumLengthValidator",
    },
    {
        "NAME": "django.contrib.auth.password_validation.CommonPasswordValidator",
    },
    {
        "NAME": "django.contrib.auth.password_validation.NumericPasswordValidator",
    },
]


DEFAULT_AUTO_FIELD = "django.db.models.BigAutoField"

USE_I18N = True

USE_L10N = True

USE_TZ = True

DATA_UPLOAD_MAX_NUMBER_FIELDS = None 

STATIC_URL = '/static/static/'
MEDIA_URL = '/static/media/'

if DEBUG:
    MEDIA_ROOT = os.path.join(BASE_DIR, 'media/')
    STATIC_ROOT = os.path.join(BASE_DIR, 'static/')
else:
    STATIC_ROOT = '/vol/web/static'
    MEDIA_ROOT = '/vol/web/media'

DEFAULT_AUTO_FIELD = 'django.db.models.BigAutoField'
