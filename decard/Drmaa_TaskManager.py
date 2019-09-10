#!/usr/bin/env python
import drmaa
import tempfile
import os
import shutil

class Drmaa_TaskManager():
    """
        This is a wrapper class to handle running a function with drmaa.
    """
    def __init__(self, tempdir=None):
        self.session = drmaa.Session()
        self.session.initialize()
        self.jobs_current = []
        self.jobs_done = []
        self.jobs_error = []
        self.tempdir = tempdir
        
    def startJob(self, command, parameters=[], num_cpus=1, partition=None, time=None ):
        """
            As it states, starts a job.
            Input:
                command: A string, with (preferably) the full path to the command to be run
                        eg '/bin/time'
                parameters (optional): A list of parameter strings
                        eg ['-n']
                num_cpus (optional): Specify how many cpus to request for this job.
                partition: (optiona): Which parition to request on call
                time (optional): How long in minutes to expect this will run. Please format as hr:min
        """
        # Create a template with our session
        job_template = self.session.createJobTemplate()
        # Set the number of CPUs (will work with SLURM only potentially)
        job_template.nativeSpecification='-c '+str(num_cpus)
        if partition:
            job_template.nativeSpecification= job_template.nativeSpecification+" -p "+partition
        if time:
            job_template.nativeSpecification= job_template.nativeSpecification+" -t "+str(time)
        # Set our command. Ideally would be given as an absolute path
        job_template.remoteCommand = command
        # Set our parameters, given in as a list
        job_template.args = parameters
        # Path where to put the stdout.
        stdout_f = tempfile.NamedTemporaryFile(dir=self.tempdir)
        # Must add a colon to the start of the path, for unclear and undocumented reasons. 
        job_template.outputPath=':'+stdout_f.name
        stderr_f = tempfile.NamedTemporaryFile(dir=self.tempdir)
        job_template.errorPath=':'+stderr_f.name
        
        working_dir = tempfile.mkdtemp(dir=self.tempdir)
        job_template.workingDir = working_dir
        
        
        job_id = self.session.runJob(job_template)
        
        self.jobs_current.append({'template': job_template,
                                  'id':     job_id,
                                  'stderr_f': stderr_f,
                                  'stdout_f': stdout_f,
                                  'working_dir': working_dir,
                                  })
    
    
    def __update__(self):
        """
            This function runs through our current jobs, drops them into done or error as needed.
        """
        done_jobs = [job for job in self.jobs_current if self.session.jobStatus(job['id'])=='done']
        self.jobs_done = self.jobs_done+done_jobs
        for job in done_jobs:
            stdout_f = job['stdout_f']
            job['output'] = stdout_f.readlines()
            stdout_f.close()

            stderr_f = job['stderr_f']
            job['error'] = stderr_f.readlines()
            stderr_f.close()
            
            shutil.rmtree(job['working_dir'])
            
            
            self.jobs_current.remove(job)
            
        error_jobs = [job for job in self.jobs_current if self.session.jobStatus(job['id'])=='failed']
        self.jobs_error = self.jobs_error + error_jobs
        for job in error_jobs:
            stdout_f = job['stdout_f']
            job['output'] = stdout_f.readlines()
            stdout_f.close()

            stderr_f = job['stderr_f']
            job['error'] = stderr_f.readlines()
            stderr_f.close()
            self.jobs_current.remove(job)
            
            
    def activeJobs(self):
        self.__update__()
        return list(self.jobs_current)
    
    def __printjob__(self,job):
        print job['template'].remoteCommand, job['template'].args, job['id'], self.session.jobStatus(job['id'])
        if 'output' in job:
            print "-- stdout --"
            for line in job['output']:
                print line
        
        if 'error' in job:
            print "-- stderr --"
            for line in job['error']:
                print line
        
    def printStatus(self):
        self.__update__()
        
        
        if len(self.jobs_done) > 0:
            print "Completed Jobs: "
            for job in self.jobs_done:
                print "_____________"
                self.__printjob__(job)
        
        if len(self.jobs_error) > 0:
            print "Jobs that failed: "
            for job in self.jobs_error:
                print "_____________"
                self.__printjob__(job)
    
        if len(self.jobs_current) > 0:
            print "Active Jobs: "
            for job in self.jobs_current:
                print "_____________"
                self.__printjob__(job)
    
        
    