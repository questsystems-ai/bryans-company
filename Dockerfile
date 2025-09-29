FROM gnina/gnina:latest

# Install tools needed for web terminal & SSH
RUN apt-get update && apt-get install -y \
    openssh-server bash curl git nano \
 && mkdir -p /var/run/sshd

# Expose port for web terminal
EXPOSE 22
EXPOSE 8888

# Keep container alive
CMD ["/usr/sbin/sshd", "-D"]
