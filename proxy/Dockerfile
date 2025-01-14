FROM nginx:stable

# Copy the Nginx configuration file into the container
COPY ./nginx.conf /etc/nginx/conf.d/default.conf

# Expose the Nginx port
EXPOSE 80

# Start Nginx
CMD ["nginx", "-g", "daemon off;"]
