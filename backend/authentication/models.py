from django.db import models

# Create your models here.
class Nitro(models.Model):
    Job = models.CharField(max_length=255)
    Sample = models.CharField(max_length=255)
    Status = models.CharField(max_length=255)
    Time = models.DateTimeField(auto_now=True)
